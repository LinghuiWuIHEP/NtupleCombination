rm -f a.out
g++ -Wall -o merge combine.cpp RecEvent.C ImpactPoints.C XCET.C `root-config --cflags` `root-config --libs`
