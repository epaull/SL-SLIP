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
package edu.umd.cs.example;

import edu.umd.cs.psl.application.inference.LazyMPEInference;
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.LazyMaxLikelihoodMPE;
import edu.umd.cs.psl.config.*
import edu.umd.cs.psl.database.DataStore
import edu.umd.cs.psl.database.Database;
import edu.umd.cs.psl.database.Partition;
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

/* 
 * The first thing we need to do is initialize a ConfigBundle and a DataStore
 */

ConfigManager cm = ConfigManager.getManager()
ConfigBundle config = cm.getBundle("basic-example")

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
def defaultPath = System.getProperty("java.io.tmpdir")
String dbpath = config.getString("dbpath", defaultPath + File.separator + "basic-example")
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

// target predicate (OPEN): does gene A influence B over a path?
m.add predicate: "influences"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// Additional CLOSED target predicates for grounded data

// data driven, expression-based from a differential analysis mutation vs not groups
m.add predicate: "diffExpr"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// data driven, expression-based, correlation in expression 
m.add predicate: "exprCorr"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// gene ontology based
m.add predicate: "goCC"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goMF"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "goBP"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// physical PPI connections: a heat-diffusion kernel distance, normalized from 0 to 1
m.add predicate: "ppi"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// test rules that predict which go category similarities are most predictive. 
m.add rule : goBP(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1
m.add rule : goCC(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1
m.add rule : goMF(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1

// test rules that predict which go category similarities are most predictive. 
m.add rule : goBP(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1
m.add rule : goCC(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1
m.add rule : goMF(A,B) & diffExpr(A,B) >> influences(A,B), weight : 1

// 'friends' also likely to be connected in network
// encode a function to make this [0,1] where directly connected things are 1
m.add rule : influences(A,B) >> distance(A,B),  weight : 3

/*
 * Finally, we define a prior on the inference predicate sl. 
 */
m.add rule: ~influences(A,B), weight: 1

/*
 * Let's see what our model looks like.
 */
println m;


/*
 * Set-up for weight learning and data input 
 */
Partition trainPart = new Partition(0);
Partition truthPart = new Partition(1);

def dir = 'data'+java.io.File.separator+'test'+java.io.File.separator;
def trainDir = dir+'train'+java.io.File.separator;

// Load static data
for (Predicate p : [gene])
{
        println "\t\t\tREADING Ground Variable " + trainDir+p.getName()+".txt";
	insert = data.getInserter(p, trainPart)
	InserterUtils.loadDelimitedData(insert, trainDir+p.getName()+".txt");
}

// load training 'truth' data. These should have a third column, 0-1 values
// 
for (Predicate p : [goCC, goMF, goBP])
{
        println "\t\t\tREADING Training Data " + trainDir+p.getName()+".txt";
	insert = data.getInserter(p, trainPart)
	InserterUtils.loadDelimitedDataTruth(insert, trainDir+p.getName()+".txt");
}
	

println "\t\t\tLoading existing sl interactions.."
insert = data.getInserter(sl, truthPart)
InserterUtils.loadDelimitedDataTruth(insert, trainDir+sl.getName()+".txt");

//////////////////////////// weight learning ///////////////////////////
println "\t\tLEARNING WEIGHTS...";

Database trainDB = data.getDatabase(trainPart, [gene, goCC, goBP, goMF] as Set);
Database truthDB = data.getDatabase(truthPart, [sl] as Set);

LazyMaxLikelihoodMPE weightLearning = new LazyMaxLikelihoodMPE(m, trainDB, truthDB, config);
weightLearning.learn();
weightLearning.close();

println "\t\tLEARNING WEIGHTS DONE";

println m

/////////////////////////// test inference //////////////////////////////////
println "\t\tINFERRING...";

def testDir = dir+'test'+java.io.File.separator;
Partition testPart = new Partition(2);
// Load static data
for (Predicate p : [gene, sl])
{
        println "\t\t\tREADING Ground Variable " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedData(insert, testDir+p.getName()+".txt");
}

// load training 'truth' data. These should have a third column, 0-1 values
// 
for (Predicate p : [goCC, goMF, goBP])
{
        println "\t\t\tREADING Training Data " + testDir+p.getName()+".txt";
	insert = data.getInserter(p, testPart)
	InserterUtils.loadDelimitedDataTruth(insert, testDir+p.getName()+".txt");
}


// don't close the sl interactions this time, but clamp everything else
Database testDB = data.getDatabase(testPart, [gene, goCC, goMF, goBP] as Set);
LazyMPEInference inference = new LazyMPEInference(m, testDB, config);
inference.mpeInference();
inference.close();

println "\t\tINFERENCE DONE";

for (GroundAtom atom : Queries.getAllAtoms(testDB, sl))
	println atom.toString() + "\t" + atom.getValue();

