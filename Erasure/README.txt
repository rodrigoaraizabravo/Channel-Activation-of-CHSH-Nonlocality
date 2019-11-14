This folder contains scripts that aid in the search of superactivation via 
tensoring of states prepared through quantum channels. These scripts are related
to the erasure channel.

This scheme is as follows: 
Alex and Bobby generate states rhoAB and rhoA'B' which they send through 
two one-qubit quantum channels (ex. E(rhoAB)=(E.Identity)(rhoAB) ). 
We ask that such channels are at least CHSH breaking so that E(rhoAB) and 
E'(rhoA'B') cannot violate the CHSH inequality. 
We then check if E(rhoAB).E'(rhoA'B') violates CHSH. 

The scrips and their descriptions are as follow: 

(1) Erasure_Utils.py: This is the main script for our protocol and implements the unidirectional 
protocols for all channels including the erasure channel.

(2) Test_Erasure.py: tests functions related to the erasure channel for the unidirectional protocol. 
Searches for superactivation as well.

(3) Bidirectional_Erasure_Utils.py: same as (1) but for the bidirectional protocol.

(4) Test_Bidirectional_Erasure.py same as (2) but for the bidirectional protocol.
