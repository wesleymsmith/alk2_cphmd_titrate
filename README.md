Using constant pH replica exchange to calculate pKa shifts
of titratable residues in ALK2.

Explore use of constant pH simulations to enhance sampling.

This project follows the workflow detailed in AMBER tutorial 18 "Constant pH MD Example
Calculating pKas for titratable side chains in HEWL" for structure setup, minimization,
heating, and equilibration.

Following equilibration, we will diverge so that we may utilize replica exchange to
allow running the production simulation over a range of pH values while allowing
the replicas to exchange structures during the production run.
To do this we first run window equilibration for each pH window. 
This is done in series starting with the two pH windows adjacent to the physiological pH window.
After this the production run is performed using replica exchange with each replica starting from
one of the previously generated windows.
This whole process, from heating through replica exchange is then repeated three more times for a total
of 4 runs, to allow additional sampling.
