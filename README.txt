This folder contains scripts that aid in the search of superactivation via 
tensoring of states prepared through quantum channels. 

This scheme is as follows: 
Alex and Bobby generate states rhoAB and rhoA'B' which they send through 
two one-qubit quantum channels (ex. E(rhoAB)=(E.Identity)(rhoAB) ). 
We ask that such channels are at least CHSH breaking so that E(rhoAB) and 
E'(rhoA'B') cannot violate the CHSH inequality. 
We then check if E(rhoAB).E'(rhoA'B') violates CHSH. 

The scrips and their descriptions are as follow: 

(1) Superactivation_Utils.py: This is the main script for our protocol and implements the bidirection 
and unidirctional protocols for all channels except the erasure channel.

(2) SingleQubit_Test.py: tests functions related to the Horodecki condition on a single qubit.

(3) Uniderectional_NumericalSearch.py: search for superactivation using the utilities on channels 
in the unidirectional protocol.

(4) Bidirectional_NumericalSearch.py same as (3) but with the bidirectional protocol.

