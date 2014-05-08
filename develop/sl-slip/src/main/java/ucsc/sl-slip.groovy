/*
 * This file is part of the Psl software.
 * Copyright 2011-2013 University of Maryland
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package src.main.ucsc;

import edu.umd.cs.psl.application.inference.MPEInference;
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.MaxLikelihoodMPE;
import edu.umd.cs.psl.application.inference.LazyMPEInference;
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.LazyMaxLikelihoodMPE;
import edu.umd.cs.psl.config.*
import edu.umd.cs.psl.database.DataStore
import edu.umd.cs.psl.database.Database;
import edu.umd.cs.psl.database.Partition;
import edu.umd.cs.psl.database.DatabasePopulator;
import edu.umd.cs.psl.database.ReadOnlyDatabase;
import edu.umd.cs.psl.database.rdbms.RDBMSDataStore
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver.Type
import edu.umd.cs.psl.model.predicate.Predicate;
import edu.umd.cs.psl.groovy.PSLModel;
import edu.umd.cs.psl.groovy.PredicateConstraint;
import edu.umd.cs.psl.groovy.SetComparison;
import edu.umd.cs.psl.model.argument.ArgumentType;
import edu.umd.cs.psl.model.argument.GroundTerm;
import edu.umd.cs.psl.model.atom.GroundAtom;
import edu.umd.cs.psl.model.function.ExternalFunction;
import edu.umd.cs.psl.ui.functions.textsimilarity.*
import edu.umd.cs.psl.ui.loading.InserterUtils;
import edu.umd.cs.psl.util.database.Queries;

import edu.umd.cs.psl.application.learning.weight.em.HardEM;

/* 
 * The first thing we need to do is initialize a ConfigBundle and a DataStore
 */

ConfigManager cm = ConfigManager.getManager()
ConfigBundle config = cm.getBundle("sl-slip")

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
def defaultPath = System.getProperty("java.io.tmpdir")
String dbpath = config.getString("dbpath", defaultPath + File.separator + "sl-slip")
println dbpath;
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, true), config)

/*
 * Now we can initialize a PslModel, which is the core component of Psl.
 * The first constructor argument is the context in which the PslModel is defined.
 * The second argument is the DataStore we will be using.
 */
PSLModel m = new PSLModel(this, data)

/* 
 * We create three predicates in the model, giving their names and list of argument types
 */
////////////////////////// predicate declaration ////////////////////////
println "\t\tDECLARING PREDICATES";

// atom predicate 
m.add predicate: "gene"        , types: [ArgumentType.UniqueID, ArgumentType.String]

// target predicate (Closed): are these genes synthetic lethal? 
m.add predicate: "slObserved"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// target predicate (OPEN): are these genes synthetic lethal? 
m.add predicate: "sl"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// Additional CLOSED target predicates for grounded data

// gene ontologies: distance scores normalized 0 to 1
m.add predicate: "goCC"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goMF"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goBP"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// physical PPI connections: a heat-diffusion kernel distance, normalized from 0 to 1
m.add predicate: "ppiConnected"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/* 
 * The 'Enemy of my Enemy is my Friend' Rule
 * This is the key network related rule: if two genes are connected by 2 hops in the sl network, they cannot be
 * sl (i.e directly connected). This could also be called the 'no triangles' rule...
 *
*/

// this doesn't seem to work: we need a set comparison: make sure these are different rules
// is ^ symmetry or not equal?
// when PSL grounds these, it won't ground the symmetrical case
// another option is to duplicate these, and add symmetry
// Experiment: reverse columns in the dataset, try both cases
m.add rule : ( sl(A,X) & sl(B,X) & (A - B) ) >> ~sl(A,B),  weight : 1


// Require a smaller intersection of the neighborhoods
// could use GO predicates for neighbors, but not for SL
// non-convex because the denominator of the Jaccard is changing. Sort of an alternating
// update rule, where we freeze the Jaccard might work. 
//m.add setcomparison: "mutuallyExclusive" , using: SetComparison.Equality, on : sl
//m.add rule :  mutuallyExclusive( {A.sl + A.sl(inv) } , {B.sl + B.sl(inv) } ) >> sl(A,B) , weight : 3

/*
 * Now we relate this to ontologies: these should be initialized by a correlation analysis of the data
 */
//m.add rule : ( sl(A,C) & sl(C,B) ) >> goCC(A,B),  weight : 1
//m.add rule : ( sl(A,C) & sl(C,B) ) >> goMF(A,B),  weight : 1
//m.add rule : ( sl(A,C) & sl(C,B) ) >> goBP(A,B),  weight : 1
m.add rule : goBP(A,B) >> sl(A,B), weight : 1
m.add rule : goCC(A,B) >> sl(A,B), weight : 1
m.add rule : goMF(A,B) >> sl(A,B), weight : 1


// 'friends' also likely to be connected in network
m.add rule : ppiConnected(A,B) >> ~sl(A,B), weight : 1


// observed values --> also SL equivalent
m.add rule : slObserved(A,B) >> sl(A,B), weight : 1

/*
 * Let's try to add an additional set comparison rule that says similarity in the negative 
 * m.add setcomparison: "connectedsl" , using: SetComparison.Equality, on : sl
 *m.add rule :  ( ~sl(A,B) >> connectedsl( {A.knows + A.knows(inv) } , {B.knows + B.knows(inv) } ) , weight : 3.2
*/

// don't need to switch data columns, this will populate the DB
m.add PredicateConstraint.Symmetric, on : sl
m.add PredicateConstraint.Symmetric, on : slObserved
m.add PredicateConstraint.Symmetric, on : goCC
m.add PredicateConstraint.Symmetric, on : goBP
m.add PredicateConstraint.Symmetric, on : goMF
m.add PredicateConstraint.Symmetric, on : ppiConnected

/*
 * Finally, we define a prior on the inference predicate sl. 
 */
m.add rule: ~sl(A,B), weight: 1

/*
 * Let's see what our model looks like.
 */
println m;


/*
 * Set-up for weight learning and data input 
 */
// everything we have that's been observed
Partition trainPart = new Partition(0);
Partition labelsPart = new Partition(1);

def dir = '../../data/test/';
def trainDir = dir+'train'+java.io.File.separator;

// Load static data
for (Predicate p : [gene, ppiConnected])
{
        println "\t\t\tREADING Ground Variable " + trainDir+p.getName()+".txt";
	insert = data.getInserter(p, trainPart)
	InserterUtils.loadDelimitedData(insert, trainDir+p.getName()+".txt");
}

// load training 'truth' data. These should have a third column, 0-1 values
// 
for (Predicate p : [slObserved, goCC, goMF, goBP])
{
        println "\t\t\tREADING Training Data " + trainDir+p.getName()+".txt";
	insert = data.getInserter(p, trainPart)
	InserterUtils.loadDelimitedDataTruth(insert, trainDir+p.getName()+".txt");
}

// observed partition
// write partition
// 	truth partition are the labels
// separate 
println "\t\t\tLoading existing sl interactions.."
insert = data.getInserter(sl, labelsPart)
InserterUtils.loadDelimitedDataTruth(insert, trainDir+sl.getName()+".txt");

//////////////////////////// weight learning ///////////////////////////
println "\t\tLEARNING WEIGHTS...";

// set up training DB -- all observed values, SL is there even though you haven't closed it
// it knows what it needs to learn, and it knows to look for pairs, but if it finds one that 
// isn't in the input file (A,X) -- populate the database method bypasses that

// PREDICATE when you have all pairs, and when you don't.
// allInteractions

// im not using labels, link prediction tasks, which is why we're running EM semi-supervised
// not a supervised task

// 
Database trainDB = data.getDatabase(trainPart, [gene, slObserved, ppiConnected, goCC, goBP, goMF] as Set);
Database labelsDB = data.getDatabase(labelsPart, [sl] as Set);

// populate database
// use this for MaxLikelihood learning (non-lazy) once it's working
DatabasePopulator dbPop = new DatabasePopulator(trainDB);
dbPop.populateFromDB(labelsDB, sl);

LazyMaxLikelihoodMPE weightLearning = new LazyMaxLikelihoodMPE(m, trainDB, labelsDB, config);
weightLearning.learn();

trainDB.close();
labelsDB.close();

//weightLearning.close();

/*
// freezing inferences as if observations
Database observedDB = data.getDatabase(readPart, [sl] as Set);

// only works because goXX has all pairs
DatabasePopulator dbPop = new DatabasePopulator(trainDB);
// this is looking at goCC, every possible ground atom it has, and take all these pairings. 
dbPop.populateFromDB(trainDB, goCC);

// now run EM on trainDB: model, random variable DB, training DB and config bundle
HardEM weightLearning = new HardEM(m, trainDB, observedDB, config);

// this might be the right way, but one equivalent would be to create a separate RV with only 
// inferred SL links

*/
println "\t\tLEARNING WEIGHTS DONE";

println m

/////////////////////////// test inference //////////////////////////////////
println "\t\tINFERRING...";

def testDir = dir+'test'+java.io.File.separator;
Partition testPart = new Partition(2);
// Load static data
for (Predicate p : [gene, ppiConnected])
{
        println "\t\t\tREADING Ground Variable " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedData(insert, testDir+p.getName()+".txt");
}

// load training 'truth' data. These should have a third column, 0-1 values
// 
for (Predicate p : [slObserved, goCC, goMF, goBP])
{
        println "\t\t\tREADING Training Data " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedDataTruth(insert, testDir+p.getName()+".txt");
}


// don't close the sl interactions this time, but clamp everything else
Database testDB = data.getDatabase(testPart, [gene, slObserved, ppiConnected, goCC, goMF, goBP] as Set);
LazyMPEInference inference = new LazyMPEInference(m, testDB, config);
inference.mpeInference();
inference.close();

println "\t\tINFERENCE DONE";

for (GroundAtom atom : Queries.getAllAtoms(testDB, sl))
	println atom.toString() + "\t" + atom.getValue();
