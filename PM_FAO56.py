# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:17:09 2019

@author: LI Baoni
"""

import numpy as np
import os,math
os.chdir('E:/YANGTZE')


#53年19358天 17年6209天 19年6940天
nsta = 59
day1 = 19358
day2 = 19358
year = 53


#注意经纬度顺序,LON first
sta = np.loadtxt('YZstation/yrd1.txt',dtype=int,skiprows=1,usecols=0).reshape(nsta) 
lat = np.loadtxt('YZstation/yrd1.txt',skiprows=1,usecols=2).reshape(nsta) 
ele = np.loadtxt('YZstation/yrd1.txt',skiprows=1,usecols=3).reshape(nsta) 
doy = np.loadtxt('YZstation/DOY1.txt',skiprows=1,usecols=1).reshape(day1)
doy8d_53 = np.loadtxt('YZstation/doy8d_53.txt',skiprows=1,usecols=1).reshape(year)
datenumber = np.loadtxt('YZstation/datenumber11.txt',usecols=1).reshape(year)


ET0 = np.zeros([nsta,day2]) 
EA = np.zeros([nsta,day2]) 
VPD = np.zeros([nsta,day2]) 
Q = np.zeros([nsta,day2]) 
RH = np.zeros([nsta,day2]) 
ET0am = np.zeros([nsta,year]) 
EAam = np.zeros([nsta,year]) 
VPDam = np.zeros([nsta,year]) 
Qam = np.zeros([nsta,year]) 
RHam = np.zeros([nsta,year]) 


for i in range(0,nsta): 
    f = str(sta[i])
    filename1 = 'E:/YRD6517/PRS/'+f+'.txt'
    filename2 = 'E:/YRD6517/RH/'+f+'.txt'
    filename3 = 'E:/YRD6517/SSD/'+f+'.txt'
    filename4 = 'E:/YRD6517/WIN/'+f+'.txt'
    filename5 = 'E:/YRD6517/TEM/'+f+'.txt'
    filename6 = 'E:/YRD6517/TEM/'+f+'.txt'
    prs = np.loadtxt(filename1,skiprows=1,usecols=2).reshape(day1)
    rh = np.loadtxt(filename2,skiprows=1,usecols=2).reshape(day1)
    ssd = np.loadtxt(filename3,skiprows=1,usecols=2).reshape(day1)
    win = np.loadtxt(filename4,skiprows=1,usecols=2).reshape(day1)
    tmax = np.loadtxt(filename5,skiprows=1,usecols=3).reshape(day1)
    tmin = np.loadtxt(filename6,skiprows=1,usecols=4).reshape(day1)
    ###########################################################################舍去2000年和2018年，需要修改
    for j in range(0,day1):
        #日地间相对距离的倒数dr
        dr = 1+0.033*math.cos(2*math.pi*doy[j]/365)
        #太阳磁偏角delta
        delta = 0.409*math.sin(2*math.pi*doy[j]/365-1.39)
        #地理纬度转换成弧度phi，北半球纬度为正
        phi = math.pi*lat[i]/180
        #日落时角ws
        ws = math.acos(-math.tan(phi)*math.tan(delta))
        #天顶辐射Ra，0.0820为太阳常数. megajules per square meter per day
        Ra = (24*60*0.0820/math.pi)*dr*(ws*math.sin(phi)*math.sin(delta)+math.cos(phi)*math.cos(delta)*math.sin(ws))
        #白昼时间N（可能最大日照时数）
        N = (24/math.pi)*ws
        #太阳辐射Rs，站点日照时数为0.1h，处理缺测值，使其等于最近日期值
        a = j
        while ssd[j]==-999 or ssd[j]==32766:
            a-= 1
            ssd[j] = ssd[a]
        Rs = (0.25+0.50*(0.1*ssd[j]/N))*Ra
        #晴空太阳辐射Rso
        Rso = (0.75+2*math.pow(10,-5)*ele[i])*Ra
        #净太阳辐射或短波辐射Rns(假设参考作物反射率为0.23)
        Rns = (1-0.23)*Rs
        #饱和水汽压es
        a = j
        while tmax[j]==-999 or tmax[j]==32766:
            a-= 1
            tmax[j] = tmax[a]
        a = j
        while tmin[j]==-999 or tmin[j]==32766:
            a-= 1
            tmin[j] = tmin[a]
        etmax = 0.6108*math.exp((17.27*0.1*tmax[j])/(0.1*tmax[j]+237.3))
        etmin = 0.6108*math.exp((17.27*0.1*tmin[j])/(0.1*tmin[j]+237.3))
        es = (etmax+etmin)/2.0
        #实际水汽压ea
        a = j
        while rh[j]==-999 or rh[j]==32766:
            a-= 1
            rh[j] = rh[a]
        ea = es*rh[j]/100
        #饱和水汽压差VPD
        a = j
        while prs[j]==-999 or prs[j]==32766:
            a-= 1
            prs[j] = prs[a]
        vpd = es-ea
        q = 0.622*ea/(prs[j]*0.01-0.378*ea)  #g/g
        #净长波辐射Rnl  中文手册里分母是4，英文是2
        Rnl = 4.903*math.pow(10,-9)*((math.pow((0.1*tmax[j]+273.16),4)+math.pow((0.1*tmin[j]+273.16),4))/4)*(0.34-0.14*math.pow(ea,0.5))*(1.35*Rs/Rso-0.35)
        #净辐射Rn
        Rn = Rns-Rnl
        #土壤热通量G
        G = 0
        #湿度计常数gamma(本站气压单位为0.1hpa，要求kpa)
        gamma = 0.665*math.pow(10,-3)*(prs[j]*0.01)
        #平均气温Tmean, 站点单位为0.1℃
        Tmean = 0.1*(tmax[j]+tmin[j])/2
        #2米风速U2，站点单位为0.1m/s
        a = j
        while win[j]==-999 or win[j]==32766:
            a-= 1
            win[j] = win[a]
        U2 = 0.1*0.748*win[j]
        #饱和水汽压曲线斜率slope
        slope = 4098*(0.6108*math.exp(17.27*Tmean/(Tmean+237.3)))/math.pow((Tmean+237.3),2)
        #参考作物腾发量ET0
        ET0[i,j] = (0.408*slope*(Rn-G)+gamma*900/(Tmean+273)*U2*(es-ea))/(slope+gamma*(1+0.34*U2)) 
        EA[i,j] = ea
        VPD[i,j] = vpd
        Q[i,j] = q
        RH[i,j] = rh[j]
        #处理计算异常值
        a = j
        while ET0[i,j]<0 or ET0[i,j]>20:
            a-= 1
            ET0[i,j] = ET0[i,a]
        a = j
        while RH[i,j]<0 or RH[i,j]>100:
            a-= 1
            RH[i,j] = RH[i,a]
        a = j
        while EA[i,j]<0 or EA[i,j]>5:
            a-= 1
            EA[i,j] = EA[i,a]
        a = j
        while VPD[i,j]<0 or VPD[i,j]>5:
            a-= 1
            VPD[i,j] = VPD[i,a]
        a = j
        while Q[i,j]<0 or Q[i,j]>5:
            a-= 1
            Q[i,j] = Q[i,a]


#折日换算成年平均
for m in range(0,nsta): #站数
    for l in range(0,year): #17年
        t = int(doy8d_53[l]) #每一年从哪一天开始
        dn = int(datenumber[l]) #每年多少天
        sum = 0
        for n in range(t,t+dn): #一年内
            sum+= ET0[m,n]
        ET0am[m,l] = sum

#for m in range(0,nsta): #站数
#    for l in range(0,year): #17年
#        t = int(doy8d_17[l]) #每一年从哪一天开始
#        dn = int(datenumber[l]) #每年多少天
#        sum = 0
#        for n in range(t,t+dn): #一年内
#            sum+= EA[m,n]
#        EAam[m,l] = sum/dn
#
#for m in range(0,nsta): #站数
#    for l in range(0,year): #17年
#        t = int(doy8d_17[l]) #每一年从哪一天开始
#        dn = int(datenumber[l]) #每年多少天
#        sum = 0
#        for n in range(t,t+dn): #一年内
#            sum+= VPD[m,n]
#        VPDam[m,l] = sum/dn
#
#for m in range(0,nsta): #站数
#    for l in range(0,year): #17年
#        t = int(doy8d_17[l]) #每一年从哪一天开始
#        dn = int(datenumber[l]) #每年多少天
#        sum = 0
#        for n in range(t,t+dn): #一年内
#            sum+= Q[m,n]
#        Qam[m,l] = sum/dn
#
#for m in range(0,nsta): #站数
#    for l in range(0,year): #17年
#        t = int(doy8d_17[l]) #每一年从哪一天开始
#        dn = int(datenumber[l]) #每年多少天
#        sum = 0
#        for n in range(t,t+dn): #一年内
#            sum+= RH[m,n]
#        RHam[m,l] = sum/dn


np.savetxt('YRD/PET.txt',ET0am,fmt="%.4f")
#np.savetxt('YRD/EA.txt',EAam,fmt="%.4f")
#np.savetxt('YRD/VPD.txt',VPDam,fmt="%.4f")
#np.savetxt('YRD/Q.txt',Qam,fmt="%.6f")
#np.savetxt('YRD/RH.txt',RHam,fmt="%.2f")


#折算成每8天
#for m in range(0,nsta): #站数
#    num = int(0)
#    for l in range(0,17): #17年
#        t = int(doy8d_17[l]) #每一年从哪一天开始
#        for k in range(0,45): #一年内从第几天开始加
#            sum = 0
#            n = t+k*8
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            n+= 1
#            sum+= ET0[m,n]
#            ET08d[m,num] = sum
#            num+= 1
#        l+= 1
#        if l==17: #最后的5或6天加一起
#            end = 6209
#        else:
#            end = int(doy8d_17[l])
#        sum = 0
#        for c in range(n+1,end):
#            sum+= ET0[m,c]
#        ET08d[m,num] = sum
#        num+= 1
#np.savetxt('AMET0/ET05.txt',ET08d,fmt="%.4f")