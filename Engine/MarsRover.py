'''
APIs for Offloading with Mars Rover
'''
from VisualizeMars import *
#from RobotCloudModel import *
from PythonRobotics.PathTracking.model_predictive_speed_and_steer_control.model_predictive_speed_and_steer_control import *
import sys,os
import numpy as np
import time
from Parameters import *
from BicycleModel import *

class Rover:
    '''
    APIs for Offloading on the Mars Rover data
    '''

    def __init__(self,env_id=ENV_ID,mars_point_cloud_fname=MARS_POINT_CLOUD):
        self.mars_point_cloud_fname=mars_point_cloud_fname
        self.env_id=env_id
        self.image=[]
        self.obs=[]
        self.loadData()

    def loadData(self):
        '''
        Load the data in objects from the file
        '''
        dataset = np.load(self.mars_point_cloud_fname)
        self.image  = dataset['image'][self.env_id]
        self.obs = dataset['point_cloud'][self.env_id]


    def getReachSets(self):
        '''
        Get Reachable Sets of Robot and Cloud
        '''
        VisualizeMarsRover.vizMarsPointClouds(self.obs,self.image)

        if self.env_id==2:
            path=[
            [(75,25),(52,70),(55,110),
            (100,100),
            (70,135),(80,190),(140,200),
            (136,160),
            (150,100),(100,50),(85,40)]
            ]
        elif self.env_id==0:
            path=[
            [(25,50),(85,70),(75,118),(100,142),(71,175),(15,167)]
            ]
        elif self.env_id==1:
            case=4
            if case==2:
                path=[
                [(20,115),(50,125),(100,129),(112,160),(165,170)]
                ]
            elif case==3:
                path=[
                [(20,115),(65,79),(100,85),(125,128),(165,170)]
                ]
            elif case==4:
                path=[
                [(20,115),(60,60),(128,76),(165,170)]
                ]
        else:
            print("Invalid Environment Chosen!")
            exit(0)
        #VisualizeMarsRover.vizMarsTraj(self.obs,self.image,path)
        (T, X, Y, yawList, V, D, A, spline_path) = doMPC(path)

        path_mpc=[]
        for x,y in zip(X,Y):
            path_mpc.append((x,y))
        VisualizeMarsRover.vizMarsTrajMPC(self.obs,self.image,path,spline_path,path_mpc)

        rsList=BicycleModel.getReachSets(X, Y, yawList, V, D, A)
        VisualizeMarsRover.vizMarsTrajRS(self.obs,self.image,spline_path,path_mpc,path,rsList)


    def getCellOrder(self):
        '''
        Get cell orders
        '''
        v=0.1
        phi=0.1
        delta=0.1
        cellOrder=BicycleModel.getCellOrderBasedError(v,phi,delta)
        return cellOrder



if True:
    rov=Rover()
    #rov.getReachSets()
    rov.getCellOrder()
