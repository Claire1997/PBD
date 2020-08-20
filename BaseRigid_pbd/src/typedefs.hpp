// ----------------------------------------------------------------------------
// typedefs.hpp
//
//  Created on: 13 Feb 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Global renaming for types
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

typedef float Real;
typedef unsigned int tIndex;

struct {
    // originally defined
    int constraintIterations = 10;
    Real deltaT = 1e-3;
    Real damping = 0.05;
    Real compressionStiffness = 1;
    Real stretchingStiffness = 1;
    Real bendingStiffness = 1;
    Real restitution = 0.01; // used for friction
    std::vector<int> staticParticles = std::vector<int>();
    
    // calculated during simulation
    int numParticles = 0;
    int numFaces = 0;
    Real originalVolume = 0;
    
    // not currently using
    Real pressure = 1;
    Real breakage = 1.5; // how much of the orig length before breaking the length
} sim;


#endif  /* _TYPEDEFS_H_ */
