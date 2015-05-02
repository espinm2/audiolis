#ifndef STATS_H
#define STATS_H

// This class was to be passed around so that we can keep track of 
// of some basic statisitics in the system.

#include <iostream>
#include <cstdio>
#include <map>
#include <vector>
#include <algorithm>

typedef unsigned int uint;

class Stats{
public:

  // Default constructor 
  Stats(){
    particle_size = 0;
  }

  // Insert function to save which shapes we would see most often
  void saveShape(uint point_num){
    shapeBin[point_num]++;
    particle_size++;
  }

  // Range in with we find these particles in
  void saveRange(double min, double max){
    minRange.push_back(min);
    maxRange.push_back(max);
  }

  // Clears the stats we are collecting
  void clearStats(){
    shapeBin.clear();
  }

  void printout(){
    
    // Sort the ranges to the medium of all of them
    std::sort(minRange.begin(), minRange.end());
    std::sort(maxRange.begin(), maxRange.end());
    int midpoint = minRange.size() / 2;
    double min = minRange[midpoint];
    double max = maxRange[midpoint];
    printf("Median Range from Min: %3.6fcm to %3.6fcm\n", min * 100, max * 100);

    // Make a print out of the map
    std::map<uint,uint>::const_iterator it;
    for( it = shapeBin.begin(); it != shapeBin.end(); it++){
      uint shape = it->first;
      uint freq  = it->second;
      double percent = freq / ((double)particle_size);
      printf("Cluster Size %d: Freq: %d ( %3.3f percent )\n",shape,freq,percent);
    }
  }

private:

  // Memeber Variables
  std::map<uint,uint>shapeBin;
  std::vector<double> minRange;
  std::vector<double> maxRange;
  int particle_size;

};


#endif