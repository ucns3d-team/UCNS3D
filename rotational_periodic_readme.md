
On the pointwise file you will set two periodic interface:
 "positive" and  "negative". The positive have an angle < respect to the negative in respect to x axis.

When you export the fluent.msh file you have to open it and search for the word "interface":
you will see two matches:

if the match is the "negative interface" change the line:
 ```
(0 "Zone 7 4868 faces 4066400..4071267, BC: negative interface = 24")
(13 (7 3e0c60 3e1f63 18 0)(
 ```
into 
 ```
(0 "Zone 7 4868 faces 4066400..4071267, BC: negative interface = 24")
(13 (7 3e0c60 3e1f63 80 0)(
 ```

so you have changed the code "18" into "80"


if the match is the "positive interface"
 ```
(0 "Zone 8 4868 faces 4071268..4076135, BC: positive interface = 24")
(13 (8 3e1f64 3e3267 18 0)(
 ```
change into 
 ```
(0 "Zone 8 4868 faces 4071268..4076135, BC: positive interface = 24")
(13 (8 3e1f64 3e3267 8 0)(
 ```
so from "18" to "8".

This was made so the code knows what is the positive and the negative interface that must be rotated of theta or -theta.

In the ROTFRAME.dat there are two new parameter:
BOUNDARY CONDITIONS: |0: Non-Periodic | 1: Periodic
1 180.0  1.0e-08

where the second are the degrees (180.0 for caradonna, 90.0 for psp and uh-60 and so on) while the second is the tolerance for the procedure that associates the periodic neighours.
you can try different value but 1.0e-08 should not give any problem. however if there are problem the simulation will crash immediately so you will know about that.