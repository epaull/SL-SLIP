
laptop-run:
	java -cp ./target/classes:`cat classpath.EvansLaptop.out` src.main.ucsc.sl-slip

toy:
	java -cp ./target/classes:`cat classpath.EvansLaptop.out` src.main.ucsc.sl-slip-freeze1

k-run:
	java -cp ./target/classes:`cat classpath.Kolossus.out` src.main.ucsc.sl-slip

k-baseline:
	java -cp ./target/classes:`cat classpath.Kolossus.out` src.main.ucsc.sl-slip.noCollab

laptop-causal:
	java -cp ./target/classes:`cat classpath.EvansLaptop.out` src.main.ucsc.causal

k-causal:
	java -cp ./target/classes:`cat classpath.Kolossus.out` src.main.ucsc.causal

laptop-compile:
	mvn compile
