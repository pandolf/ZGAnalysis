#! /usr/bin/python

if __name__ == '__main__':

  import os
  import sys
  from optparse import OptionParser
  parser = OptionParser()
  parser.usage = ""
  parser.add_option("-s","--significance", dest="significance",
                    default="",
                    help="compute significance instead of limit")
  parser.add_option("-t","--toys", dest="toys",
                    default="",
                    help="use toys instead of asymptotic")


  (options,args) = parser.parse_args()

  thisdir = args[0]
  while thisdir.endswith("/"): thisdir = thisdir[:-1]

  lines = [] 

  for d in os.listdir(thisdir) :

    if not d.startswith("datacard"): continue
    signalname = d.split(thisdir)[1]
    signalname = signalname.split(".txt")[0]
    signalname = signalname[1:]

    print "Starting: " + signalname
    mass = signalname.split("test")[1]

    logfilename = thisdir + "/limit_" + signalname + ".log"
    os.system("combine -M Asymptotic -L libdiphotonsUtils -m " + str(mass) + " -d " + thisdir + "/" + d + " -t -1 >& " + logfilename )

    logfile = open( logfilename, "r" )

    obs = 0.
    exp = 0.
    exp_p1s = 0.
    exp_p2s = 0.
    exp_m1s = 0.
    exp_m2s = 0.

    for line in logfile:

      if not (line.startswith("Expected") or line.startswith("Observed")) : continue

      limit = line.split("r < ")[1]
      limit.rstrip()
      if "Observed"       in line: obs     = float(limit)
      if "Expected  2.5%" in line: exp_m2s = float(limit)
      if "Expected 16.0%" in line: exp_m1s = float(limit)
      if "Expected 50.0%" in line: exp     = float(limit)
      if "Expected 84.0%" in line: exp_p1s = float(limit)
      if "Expected 97.5%" in line: exp_p2s = float(limit)

    lines.append( "m: " + str(mass) + " obs: " + str(obs) + " exp: " + str(exp) + " exp_m1s: " + str(exp_m1s) + " exp_m2s: " + str(exp_m2s) + " exp_p1s: " + str(exp_p1s) + " exp_p2s: " + str(exp_p2s) + "\n" )

  lines.sort()

  outfile = open( thisdir + "/limits.txt", "w" )
  for line in lines: outfile.write(line)

  
