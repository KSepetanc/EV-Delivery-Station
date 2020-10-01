* Txt script
$onecho > input_data_copy.txt
i=data.xlsx
o=data.gdx

par=N_start              rng=N_start!a2                  Rdim=1 Cdim=0
par=Eta                  rng=Eta!a1:v22                  Rdim=1 Cdim=1
par=Lambda_week          rng=Price!b2                    Rdim=2 Cdim=0
par=N_out                rng=N_out!b2                    Rdim=3 Cdim=0
$offecho

* Create gdx from excel using txt script
$call gdxxrw @input_data_copy.txt

* Delete txt script
$call del /f input_data_copy.txt




















