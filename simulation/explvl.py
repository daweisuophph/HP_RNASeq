import re
import os
import math
import numpy.random as rnd
from sets import Set

numOfSubjects = 5

mu = -2
sigma = 2

test_folder_name = "reads_repeated_exp"
data_path = "./"
alpha = 1

for exp_id in xrange(1, 11):
   rnd.seed(0) # reset the seed so the explvl is the same
   test_name = test_folder_name + "/exp"+ str(exp_id);
   for s in range(1,numOfSubjects+1):
      name = data_path+test_name + "/"+str(s)+"N/"
      d = os.path.dirname(name)
      if not os.path.exists(d):
         os.makedirs(d)
      name = data_path+test_name + "/"+str(s)+"T/"
      d = os.path.dirname(name)
      if not os.path.exists(d):
         os.makedirs(d)

   # NON_DE
   fp = [None] * (numOfSubjects+1)
   for s in range(1, numOfSubjects+1):
      fp[s] = [None]*2
      fp[s][0] = open(data_path+test_name+"/"+str(s)+"N/explvprofileNONDE.txt", "w")
      fp[s][1] = open(data_path+test_name+"/"+str(s)+"T/explvprofileNONDE.txt", "w")

   for s in range(1, numOfSubjects+1):
      fp[s][0].write('#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup\tExplv\n')
      fp[s][1].write('#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup\tExplv\n')

   f = open('NON_DE/template.txt', 'r')
   line = f.readline()
   line = f.readline()
   og = -1
   iso = 0
   while line != '':
      iso = iso + 1
      features = re.split('[\t\n]+', line)
      g = int(features[5])
      if g != og:
         #print g
         og = g
         K = int(features[6])
         iso = 0
         
         n = rnd.lognormal(mu, sigma)
         #n = 1000
         p = rnd.dirichlet((alpha,)*K, 1)
         p = p * n

      levelT = p[0][iso] * int(features[1]);
      levelN = p[0][iso] * int(features[1]);

      #print str(levelT) + " " + str(levelN)

      for s in range(1, numOfSubjects+1):
         lvl = rnd.normal(levelN, levelN*0.1)
         while lvl < 0:
            lvl = rnd.normal(levelN, levelN*0.1)
         fp[s][0].write(features[0]+"\t"+features[1]+"\t"+features[2]+"\t"+features[3]+"\t" \
               +features[4]+"\t"+features[5]+"\t"+features[6]+"\t"+str(lvl)+"\n")

         lvl = rnd.normal(levelT, levelT*0.1)
         while lvl < 0:
            lvl = rnd.normal(levelT, levelT*0.1)
         fp[s][1].write(features[0]+"\t"+features[1]+"\t"+features[2]+"\t"+features[3]+"\t" \
               +features[4]+"\t"+features[5]+"\t"+features[6]+"\t"+str(lvl)+"\n")
      
      line = f.readline()
   f.close()

   for s in range(1, numOfSubjects+1):
      fp[s][0].close()
      fp[s][1].close()

   # DE
   fp = [None] * (numOfSubjects+1)
   for s in range(1, numOfSubjects+1):
      fp[s] = [None]*2
      fp[s][0] = open(data_path+test_name+"/"+str(s)+"N/explvprofileDE.txt", "w")
      fp[s][1] = open(data_path+test_name+"/"+str(s)+"T/explvprofileDE.txt", "w")

   for s in range(1, numOfSubjects+1):
      fp[s][0].write('#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup\tExplv\n')
      fp[s][1].write('#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup\tExplv\n')

   f = open('DE/template.txt', 'r')
   line = f.readline()
   line = f.readline()
   og = -1
   iso = 0
   while line != '':
      iso = iso + 1
      features = re.split('[\t\n]+', line)
      g = int(features[5])
      if g != og:
         #print g
         og = g
         K = int(features[6])
         iso = 0
         
         n = rnd.lognormal(mu, sigma)
         #n = 1000
         p = rnd.dirichlet((alpha,)*K, 1)
         q = rnd.dirichlet((alpha,)*K, 1)
         p = p * n
         q = q * n

      levelN = p[0][iso] * int(features[1]);
      levelT = q[0][iso] * int(features[1]);

      #print str(levelT) + " " + str(levelN)

      for s in range(1, numOfSubjects+1):
         lvl = rnd.normal(levelN, levelN*0.1)
         while lvl < 0:
            lvl = rnd.normal(levelN, levelN*0.1)
         fp[s][0].write(features[0]+"\t"+features[1]+"\t"+features[2]+"\t"+features[3]+"\t" \
               +features[4]+"\t"+features[5]+"\t"+features[6]+"\t"+str(lvl)+"\n")

         lvl = rnd.normal(levelT, levelT*0.1)
         while lvl < 0:
            lvl = rnd.normal(levelT, levelT*0.1)
         fp[s][1].write(features[0]+"\t"+features[1]+"\t"+features[2]+"\t"+features[3]+"\t" \
               +features[4]+"\t"+features[5]+"\t"+features[6]+"\t"+str(lvl)+"\n")
      
      line = f.readline()
   f.close()

   for s in range(1, numOfSubjects+1):
      fp[s][0].close()
      fp[s][1].close()
