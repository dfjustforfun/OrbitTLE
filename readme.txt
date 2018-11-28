
/etc/profile

# set OrbitTLE environment
export PYTHONPATH="/home/cpf/Documents/OrbitTLE:$PYTHONPATH"


e.g.
python3 calculateOrbitTLE.py --lon 113.7783 --lat 34.7444 --alt 0.07 --startTime 2018-08-21/0:0:0 --endTime 2018-08-21/23:59:59

TODO:
1. The SDP4 process needs validation.

### 2018-11-28

1. change the format of output data
