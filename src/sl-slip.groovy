/*
 * This file is part of the PSL software.
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

/*
 * A ConfigBundle is a set of key-value pairs containing configuration options. One place these
 * can be defined is in psl-example/src/main/resources/psl.properties
 */
ConfigManager cm = ConfigManager.getManager()
ConfigBundle config = cm.getBundle("basic-example")

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
def defaultPath = System.getProperty("java.io.tmpdir")
String dbpath = config.getString("dbpath", defaultPath + File.separator + "basic-example")
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, true), config)

/*
 * Now we can initialize a PSLModel, which is the core component of PSL.
 * The first constructor argument is the context in which the PSLModel is defined.
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


// target predicate (OPEN): are these genes synthetic lethal? 
m.add predicate: "SL"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]


// Additional CLOSED target predicates for grounded data

// gene ontologies
m.add predicate: "go_CC"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "go_MP"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "go_BP"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// physical PPI connections
m.add predicate: "ppi"	, types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/* 
 * The 'Enemy of my Enemy is my Friend' Rule
 * This is the key network related rule: if two genes are connected by 2 hops in the SL network, they cannot be
 * SL (i.e directly connected). This could also be called the 'no triangles' rule...
 *
*/
m.add rule : ( SL(A,X) & SL(B,Y) ) >> ~SL(A,B),  weight : 5

/*
 * Now we relate this to ontologies: these should be initialized by a correlation analysis of the data
 */
m.add rule : ( SL(A,X) & SL(B,Y) ) >> go_CC(A,B),  weight : 1
m.add rule : ( SL(A,X) & SL(B,Y) ) >> go_MF(A,B),  weight : 1
m.add rule : ( SL(A,X) & SL(B,Y) ) >> go_BP(A,B),  weight : 1

// 'friends' also likely to be connected in network
m.add rule : ( SL(A,X) & SL(B,Y) ) >> ppi(A,B),  weight : 3

/*
 * Let's try to add an additional set comparison rule that says similarity in the negative 
 * m.add setcomparison: "connectedSL" , using: SetComparison.Equality, on : SL
 *m.add rule :  ( ~SL(A,B) >> connectedSL( {A.knows + A.knows(inv) } , {B.knows + B.knows(inv) } ) , weight : 3.2
*/

/* Next, we define some constraints for our model. In this case, we restrict that each person can be aligned to at most one other person
 * in the other social network. To do so, we define two partial functional constraints where the latter is on the inverse.
 * We also say that samePerson must be symmetric, i.e., samePerson(p1, p2) == samePerson(p2, p1).
 */
m.add PredicateConstraint.Symmetric, on : SL

/*
 * Finally, we define a prior on the inference predicate SL. 
 */
m.add rule: ~SL(A,B), weight: 1

/*
 * Let's see what our model looks like.
 */
println m;


