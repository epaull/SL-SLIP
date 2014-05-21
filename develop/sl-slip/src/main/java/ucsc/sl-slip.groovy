/*
 * SL-SLIP (Soft Logic - Synthetic Lethal Interacting Pairs)
 * This software package is designed to predict gene pairs that have 
 * synthetic lethal (SL) or synthetic sick (SS) interactions. It is built 
 * on top of the Proablistic Soft Logic (PSL) framework that allows for collective
 * inference of the joint graph-identification problem posed by this task. 
 *   
 * Template code from the PSL project was used to write this software: the
 * copyright and license for that code is below.
 *
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

// The single atom/node for this project
m.add predicate: "gene"        , types: [ArgumentType.UniqueID, ArgumentType.String]

// target predicate (Closed): are these genes synthetic lethal? 
m.add predicate: "slObserved"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// target predicate (OPEN): are these genes synthetic lethal? 
m.add predicate: "sl"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// dummy variable for populating the DB: this must include the same set of pairs as 'sl'
m.add predicate: "consider"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// Additional CLOSED target predicates for grounded data:

// gene ontologies: distance scores normalized 0 to 1
m.add predicate: "goCC"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goMF"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goBP"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// direct connections in the PPI net
m.add predicate: "ppiEdges"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// KERNEL arguments
m.add predicate: "ppiKernel"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
// G- SL kernel from Qi/Bader 2008
m.add predicate: "negKernel"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/* 
 * The 'Enemy of my Enemy is my Friend' Rule
 * This is the key network related rule: if two genes are connected by 2 hops in the sl network, they cannot be
 * sl (i.e directly connected). This could also be called the 'no triangles' rule...
 * 
*/
m.add rule : ( consider(A,B) & sl(A,X) & sl(X,B) & (A - B) ) >> ~sl(A,B),  weight : 10
// this is the 2-hop SL-PPI rule
m.add rule : ( consider(A,B) & sl(A,X) & ppiEdges(X,B) & (A - B) ) >> sl(A,B),  weight : 10

// simple hypotheses for Gene Ontology interactions: we might need to update the prior weights here
m.add rule : ( consider(A,B) & goBP(A,B) & (A-B) ) >> sl(A,B), weight : 1
m.add rule : ( consider(A,B) & goCC(A,B) & (A-B) ) >> sl(A,B), weight : 1
m.add rule : ( consider(A,B) & goMF(A,B) & (A-B) ) >> ~sl(A,B), weight : 1

// SL means generally not connected in the PPI net
m.add rule : ( consider(A,B) & ppiKernel(A,B) ) >> ~sl(A,B), weight : 10
// this should be predictive of SL, according to the Qi/Bader 2008 paper
m.add rule : ( consider(A,B) & negKernel(A,B) ) >> sl(A,B), weight : 10

// 'friends' also likely to be connected in network
m.add rule : ( consider(A,B) & ppiEdges(A,B) ) >> ~sl(A,B), weight : 10

// observed values --> also SL equivalent: this weight needs to be high
// enough to 'lock down' the sl predicate when an edge is observed
m.add rule : ( consider(A, B) & slObserved(A,B) ) >> sl(A,B), weight : 100

// all the predicates are symmetric!
m.add PredicateConstraint.Symmetric, on : sl
m.add PredicateConstraint.Symmetric, on : slObserved
m.add PredicateConstraint.Symmetric, on : goCC
m.add PredicateConstraint.Symmetric, on : goBP
m.add PredicateConstraint.Symmetric, on : goMF
m.add PredicateConstraint.Symmetric, on : ppiEdges
m.add PredicateConstraint.Symmetric, on : ppiKernel
m.add PredicateConstraint.Symmetric, on : negKernel
m.add PredicateConstraint.Symmetric, on : consider

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

def dir = '../../data/yeast/TEST/';
def trainDir = dir+'train'+java.io.File.separator;

// Load static data
for (Predicate p : [gene, consider])
{
        println "\t\t\tREADING Ground Variable " + trainDir+p.getName()+".txt";
	insert = data.getInserter(p, trainPart)
	InserterUtils.loadDelimitedData(insert, trainDir+p.getName()+".txt");
}

// load training 'truth' data. These should have a third column, 0-1 values
// 
for (Predicate p : [slObserved, goCC, goMF, goBP, ppiEdges, ppiKernel, negKernel])
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
Database trainDB = data.getDatabase(trainPart, [gene, consider, slObserved, ppiEdges, goCC, goBP, goMF, ppiKernel, negKernel] as Set);
Database labelsDB = data.getDatabase(labelsPart, [sl] as Set);

// populate database
// use this for MaxLikelihood learning (non-lazy) once it's working
DatabasePopulator dbPop = new DatabasePopulator(trainDB);
dbPop.populateFromDB(labelsDB, sl);

MaxLikelihoodMPE weightLearning = new MaxLikelihoodMPE(m, trainDB, labelsDB, config);
//HardEM weightLearning = new HardEM(m, trainDB, labelsDB, config);
weightLearning.learn();

trainDB.close();
labelsDB.close();
//weightLearning.close();

println "\t\tLEARNING WEIGHTS DONE";

println m

/////////////////////////// test inference //////////////////////////////////
println "\t\tINFERRING...";

def testDir = dir+'test'+java.io.File.separator;
Partition testPart = new Partition(2);
// Load static data without edge weights
for (Predicate p : [gene, consider])
{
        println "\t\t\tREADING Ground Variable " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedData(insert, testDir+p.getName()+".txt");
}

// load training 'truth' data. These are closed predicates with a third column, 0-1 values
// 
for (Predicate p : [slObserved, goCC, goMF, goBP, ppiEdges])
{
        println "\t\t\tREADING Training Data " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedDataTruth(insert, testDir+p.getName()+".txt");
}

// don't close the sl interactions this time, but clamp everything else
Database testDB = data.getDatabase(testPart, [gene, consider, slObserved, ppiEdges, goCC, goMF, goBP, ppiKernel, negKernel] as Set);
MPEInference inference = new MPEInference(m, testDB, config);

//create dummy partition
Partition dummyPart = new Partition(3);
insert = data.getInserter(sl, dummyPart)
InserterUtils.loadDelimitedDataTruth(insert, testDir+sl.getName()+".txt");
Database dummyDB = data.getDatabase(dummyPart, [sl] as Set);
// initialize the random variables that will get inferred with 'sl' from the learning step
DatabasePopulator dbPop2 = new DatabasePopulator(testDB);
dbPop2.populateFromDB(dummyDB, sl);

inference.mpeInference();
inference.close();

println "\t\tINFERENCE DONE";

for (GroundAtom atom : Queries.getAllAtoms(testDB, sl))
	println atom.toString() + "\t" + atom.getValue();
