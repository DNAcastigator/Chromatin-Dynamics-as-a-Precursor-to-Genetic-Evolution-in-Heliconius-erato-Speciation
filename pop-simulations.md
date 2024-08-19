# SLiM simulations
This is a general template used to run each simulation in SLiM (https://messerlab.org/slim/), for simulating Fst behavior in diverging populations, for information about each function used go to the official Manual:

We simulated two populations, each with two chromosomes. The first chromosome comprised 2% of nucleotides under weak selective forces (s=0.01),
with the remainder being neutral (s=0). The second chromosome was entirely neutral. We used this second chromosome as a control because the neutral nucleotides on the first chromosome might be subject to hitchhiking, whereas those on the second chromosome should not.
The simulations were scaled down by a factor of 100 for practical reasons, following the approach explained on page 145 of the SLiM manual.

The first step is to scale the simulation and initialize all the elements.
```eidos
//scale 100 time
initialize() {
    // Initialize mutation rates and types
    initializeMutationRate(2.9e-7); // scaled mutation rate
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", 0.5, "f", 0.01);

    // Define genomic element types
    initializeGenomicElementType("g1", c(m1,m2), c(0.98,0.02));
    initializeGenomicElementType("g2", m1, 1);

    // Define genomic elements
    initializeGenomicElement(g1, 0, 10000);
    initializeGenomicElement(g2, 10500, 20000);

    // Initialize recombination rate
    initializeRecombinationRate(2E-6);

    // Initialize tree-sequence recording

}
```

Then we can start the simulation by creating two populations with a very high gene flow
```eidos
1 early() {
    // Create populations
    sim.addSubpop("p1", 10000);
    sim.addSubpopSplit("p2", 10000,p1);
    // Set migration rates
    p1.setMigrationRates(p2, 0.5);
    p2.setMigrationRates(p1, 0.5);

}
```
we let the simulation run for 1k generations and then started to record the Fst for the different elements:

```eidos
1000 early() {

    // Create a log file and add logging cycles
    log = community.createLogFile("sim_log_mig.txt", logInterval=100);
    log.addCycle();
   log.addCustomColumn("FST_sel", "calcFST(p1.genomes, p2.genomes, c(p1.genomes.mutationsOfType(m2),p2.genomes.mutationsOfType(m2)));");
   log.addCustomColumn("FST_back", "calcFST(p1.genomes, p2.genomes, c(p1.genomes.mutationsOfType(m1),p2.genomes.mutationsOfType(m1)),0,10000);");
   log.addCustomColumn("FST_neutr", "calcFST(p1.genomes, p2.genomes, NULL,10500, 20000);");
}
```

after 2k generations we reduced the gene flow (in this example to 0) and we let the generation go for a total of 30k generations.
```eidos
2000 early() {
    // Stop migration
    p1.setMigrationRates(p2, 0);
    p2.setMigrationRates(p1, 0);


}



30000 late() {

}
```
