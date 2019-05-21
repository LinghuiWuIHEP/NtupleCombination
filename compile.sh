rm -f a.out
g++ -Wall combine.cpp RecEvent.C AhTree.C ImpactPoints.C XCET.C `root-config --cflags` `root-config --libs`
