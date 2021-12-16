# All the Parameters go here.

'''
Parameters required for the code
'''
import os,sys

'''
Please add the following line in ~/.bashrc
export ULS_ROOT_DIR = <YOUR PROJECT ROOT>
'''

PROJECT_ROOT = os.environ['ULS_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

LIB_PATH=PROJECT_ROOT+'/'+'lib/'
SRC_PATH=PROJECT_ROOT+'/'+'src/'
OUTPUT_PATH=PROJECT_ROOT+'/'+'output/'
PICKLE_PATH=PROJECT_ROOT+'/'+'pickles/'
DATA_PATH=PROJECT_ROOT+'/'+'data/'

ROVER_RESULTS=OUTPUT_PATH+'/mars_rover_results/'
MARS_POINT_CLOUD=PROJECT_ROOT+'/data/mars_data/mars_point_cloud.npz'
PER_COVERAGE=100
ENV_ID=0
RS_STEPS=1

PKPD_RESULTS=OUTPUT_PATH+'/pkpd_results/'


# Detailed Parameters
BIGM=1e500
EPSILON=1e-3
PRED_EP=1e-50

PRINT_COUNT=5
P=1

CLOSE=0.2
PRECISION=1e15

h=0.01

PRECISION=1e55

NO_SAMPPLES=50

PRED_EP=1e-3
INTERVAL=1000000000000000000
RED_INT_INTRVL=4
SAMPLES=20
VLENGTH=17

INTERVAL=99
