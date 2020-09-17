////////////////////////////////////////////////////////////////////////////////
//
//  Source code for the paper
//
//  "Fast Ray--Tetrahedron Intersection"
//
//  (c) Nikos Platis, Theoharis Theoharis 2003
//
//  Department of Informatics and Telecommunications,
//  University of Athens, Greece
//  {nplatis|theotheo}@di.uoa.gr
//
////////////////////////////////////////////////////////////////////////////////
#include <lm/core.h>
#include <lm/lm.h>

#include <glm/glm.hpp>
#include <lm/math.h>


// Pluecker ////////////////////////////////////////////////////////////////////

class Pluecker
{
private:
    lm::Vec3 d;   // direction
    lm::Vec3 c;   // cross

public:
    Pluecker()
    { };
    
    // Pluecker coordinates of a ray specified by points orig and dest.
    Pluecker(const lm::Vec3& orig, const lm::Vec3& dest)
        : d(dest-orig), c(glm::cross(dest,orig))
    { };

    // Pluecker coordinates of a ray specified by point orig and
    // direction dir. The bool parameter serves only to discriminate this
    // from the previous constructor
    Pluecker(const lm::Vec3& orig, const lm::Vec3& dir, bool haveDir)
        : d(dir), c(glm::cross(dir,orig))
    { };

    // Permuted inner product operator
    friend double operator*(const Pluecker& pl1, const Pluecker& pl2);
};


inline double operator*(const Pluecker& pl1, const Pluecker& pl2)
{
    return (glm::dot(pl1.d , pl2.c) + glm::dot(pl2.d , pl1.c));
}


////////////////////////////////////////////////////////////////////////////////
//
//  Source code for the paper
//
//  "Fast Ray--Tetrahedron Intersection"
//
//  (c) Nikos Platis, Theoharis Theoharis 2003
//
//  Department of Informatics and Telecommunications,
//  University of Athens, Greece
//  {nplatis|theotheo}@di.uoa.gr
//
////////////////////////////////////////////////////////////////////////////////

// Computes the parametric distance tEnter and tLeave 
// of enterPoint and leavePoint from orig
// (To enhance performance of the intersection algorithm, this code
// should in practice be incorporated to the function below.)

void ComputeParametricDist(
    const lm::Vec3& orig, const lm::Vec3& dir,
    const lm::Vec3& enterPoint, const lm::Vec3& leavePoint,
    double& tEnter, double& tLeave)
{
    if (dir.x)  {
        double invDirx = 1.0 / dir.x;
        tEnter = (enterPoint.x - orig.x) * invDirx;
        tLeave = (leavePoint.x - orig.x) * invDirx;
    }
    else if (dir.y)  {
        double invDiry = 1.0 / dir.y;
        tEnter = (enterPoint.y - orig.y) * invDiry;
        tLeave = (leavePoint.y - orig.y) * invDiry;
    }
    else  {
        double invDirz = 1.0 / dir.z;
        tEnter = (enterPoint.z - orig.z) * invDirz;
        tLeave = (leavePoint.z - orig.z) * invDirz;
    }
}



// Ray--tetrahedron intersection algorithm using Pluecker coordinates

bool RayTetraPluecker(
    const lm::Vec3& orig, const lm::Vec3& dir,
    const lm::Vec3 vert[],
    int& enterFace, int& leaveFace,
    lm::Vec3& enterPoint, lm::Vec3& leavePoint,
    double& uEnter1, double& uEnter2, double& uLeave1, double& uLeave2,
    double& tEnter, double& tLeave)
{
    enterFace = -1;
    leaveFace = -1;

    double uAB = 0, uAC = 0, uDB = 0, uDC = 0, uBC = 0, uAD = 0;
                    // Keep the compiler happy about uninitialized variables
    int signAB = -2, signAC = -2, signDB = -2, 
        signDC = -2, signBC = -2, signAD = -2;

    // In the following: A,B,C,D=vert[i], i=0,1,2,3.
    Pluecker plRay(orig, dir, true);

    int nextSign = 0;

    // Examine face ABC
    uAB = plRay * Pluecker(vert[0], vert[1]);
    signAB = glm::sign(uAB);

    uAC = plRay * Pluecker(vert[0], vert[2]);
    signAC = glm::sign(uAC);
    
    if ((signAC == -signAB)  ||  (signAC == 0)  ||  (signAB == 0))  {
        // Face ABC may intersect with the ray
        uBC = plRay * Pluecker(vert[1], vert[2]);
        signBC = glm::sign(uBC);

        int signABC = signAB;
        if (signABC == 0)  {
            signABC = -signAC;
            if (signABC == 0)  {
                signABC = signBC;
            }
        }

        if ((signABC != 0)  && 
            ((signBC == signABC)  ||  (signBC == 0)))  {
            // Face ABC intersects with the ray
            double invVolABC = 1.0 / (uAB + uBC - uAC);
            if (signABC == 1)  {
                enterFace = 3;
                uEnter1 = -uAC * invVolABC;
                uEnter2 =  uAB * invVolABC;
                enterPoint =   (1-uEnter1-uEnter2)*vert[0] 
                             + uEnter1*vert[1] + uEnter2*vert[2];
                
                nextSign = -1;
            }
            else  {
                leaveFace = 3;
                uLeave1 = -uAC * invVolABC;
                uLeave2 =  uAB * invVolABC;
                leavePoint =   (1-uLeave1-uLeave2)*vert[0] 
                             + uLeave1*vert[1] + uLeave2*vert[2];

                nextSign = 1;
            }

            // Determine the other intersecting face between BAD, CDA, DCB
            // Examine face BAD
            uAD = plRay * Pluecker(vert[0], vert[3]);
            signAD = glm::sign(uAD);

            if ((signAD == nextSign)  ||  (signAD == 0))  {
                // Face BAD may intersect with the ray
                uDB = plRay * Pluecker(vert[3], vert[1]);
                signDB = glm::sign(uDB);

                if ((signDB == nextSign)  ||  
                    ((signDB == 0)  &&  
                     ((signAD != 0)  ||  (signAB != 0))))  {
                    // Face BAD intersects with the ray
                    double invVolBAD = 1.0 / (uAD + uDB - uAB);
                    if (nextSign == 1)  {
                        enterFace = 2;
                        uEnter1 =  uDB * invVolBAD;
                        uEnter2 = -uAB * invVolBAD;
                        enterPoint =   (1-uEnter1-uEnter2)*vert[1]
                                     + uEnter1*vert[0] + uEnter2*vert[3];
                        
                        ComputeParametricDist(orig, dir,
                                              enterPoint, leavePoint,
                                              tEnter, tLeave);
                        return true;
                    }
                    else  {
                        leaveFace = 2;
                        uLeave1 =  uDB * invVolBAD;
                        uLeave2 = -uAB * invVolBAD;
                        leavePoint =   (1-uLeave1-uLeave2)*vert[1]
                                     + uLeave1*vert[0] + uLeave2*vert[3];

                        ComputeParametricDist(orig, dir,
                                              enterPoint, leavePoint,
                                              tEnter, tLeave);
                        return true;
                    }
                }           
            }

            // Face BAD does not intersect with the ray.
            // Determine the other intersecting face between CDA, DCB
            uDC = plRay * Pluecker(vert[3], vert[2]);
            signDC = glm::sign(uDC);

            if ((signDC == -nextSign)  ||
                ((signDC == 0)  &&  ((signAD != 0)  ||  (signAC != 0))))  {
                // Face CDA intersects with the ray
                double invVolCDA = 1.0 / (uAC - uDC - uAD);
                if (nextSign == 1)  {
                    enterFace = 1;
                    uEnter1 =  uAC * invVolCDA;
                    uEnter2 = -uDC * invVolCDA;
                    enterPoint =   (1-uEnter1-uEnter2)*vert[2] 
                                 + uEnter1*vert[3] + uEnter2*vert[0];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
                else  {
                    leaveFace = 1;
                    uLeave1 =  uAC * invVolCDA;
                    uLeave2 = -uDC * invVolCDA;
                    leavePoint =   (1-uLeave1-uLeave2)*vert[2] 
                                 + uLeave1*vert[3] + uLeave2*vert[0];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
            }
            else  {
                // Face DCB intersects with the ray
                if (signDB == -2)  {
                    uDB = plRay * Pluecker(vert[3], vert[1]);
                }

                double invVolDCB = 1.0 / (uDC - uBC - uDB);
                if (nextSign == 1)  {
                    enterFace = 0;
                    uEnter1 = -uDB * invVolDCB;
                    uEnter2 =  uDC * invVolDCB;
                    enterPoint =   (1-uEnter1-uEnter2)*vert[3] 
                                 + uEnter1*vert[2] + uEnter2*vert[1];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
                else  {
                    leaveFace = 0;
                    uLeave1 = -uDB * invVolDCB;
                    uLeave2 =  uDC * invVolDCB;
                    leavePoint =   (1-uLeave1-uLeave2)*vert[3] 
                                 + uLeave1*vert[2] + uLeave2*vert[1];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
            }
        }
    }

    // Examine face BAD
    uAD = plRay * Pluecker(vert[0], vert[3]);
    signAD = glm::sign(uAD);

    if ((signAD == -signAB)  ||  (signAB == 0)  ||  (signAD == 0))  {
        // Face BAD may intersect with the ray
        uDB = plRay * Pluecker(vert[3], vert[1]);
        signDB = glm::sign(uDB);

        int signBAD = -signAB;
        if (signBAD == 0)  {
            signBAD = signAD;
            if (signBAD == 0)  {
                signBAD = signDB;
            }
        }

        if ((signBAD != 0)  &&
            ((signDB == signBAD)  ||  (signDB == 0)))  {
            // Face BAD intersects with the ray
            double invVolBAD = 1.0 / (uAD + uDB - uAB);
            if (signBAD == 1)  {
                enterFace = 2;
                uEnter1 =  uDB * invVolBAD;
                uEnter2 = -uAB * invVolBAD;
                enterPoint =   (1-uEnter1-uEnter2)*vert[1]
                             + uEnter1*vert[0] + uEnter2*vert[3];

                nextSign = -1;
            }
            else  {
                leaveFace = 2;
                uLeave1 =  uDB * invVolBAD;
                uLeave2 = -uAB * invVolBAD;
                leavePoint =   (1-uLeave1-uLeave2)*vert[1]
                             + uLeave1*vert[0] + uLeave2*vert[3];

                nextSign = 1;
            }

            // Determine the other intersecting face between CDA, DCB
            uDC = plRay * Pluecker(vert[3], vert[2]);
            signDC = glm::sign(uDC);

            if ((signDC == -nextSign)  ||
                ((signDC == 0)  &&  ((signAD != 0)  ||  (signAC != 0))))  {
                // Face CDA intersects with the ray
                double invVolCDA = 1.0 / (uAC - uDC - uAD);
                if (nextSign == 1)  {
                    enterFace = 1;
                    uEnter1 =  uAC * invVolCDA;
                    uEnter2 = -uDC * invVolCDA;
                    enterPoint =   (1-uEnter1-uEnter2)*vert[2] 
                                 + uEnter1*vert[3] + uEnter2*vert[0];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
                else  {
                    leaveFace = 1;
                    uLeave1 =  uAC * invVolCDA;
                    uLeave2 = -uDC * invVolCDA;
                    leavePoint =   (1-uLeave1-uLeave2)*vert[2] 
                                 + uLeave1*vert[3] + uLeave2*vert[0];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
            }
            else  {
                // Face DCB intersects with the ray
                if (signBC == -2)  {
                    uBC = plRay * Pluecker(vert[1], vert[2]);
                }

                double invVolDCB = 1.0 / (uDC - uBC - uDB);
                if (nextSign == 1)  {
                    enterFace = 0;
                    uEnter1 = -uDB * invVolDCB;
                    uEnter2 =  uDC * invVolDCB;
                    enterPoint =   (1-uEnter1-uEnter2)*vert[3] 
                                 + uEnter1*vert[2] + uEnter2*vert[1];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
                else  {
                    leaveFace = 0;
                    uLeave1 = -uDB * invVolDCB;
                    uLeave2 =  uDC * invVolDCB;
                    leavePoint =   (1-uLeave1-uLeave2)*vert[3] 
                                 + uLeave1*vert[2] + uLeave2*vert[1];

                    ComputeParametricDist(orig, dir,
                                          enterPoint, leavePoint,
                                          tEnter, tLeave);
                    return true;
                }
            }
        }
    }

    // Examine face CDA
    if ((-signAD == signAC)  ||  (signAC == 0)  ||  (signAD == 0))  {
        // Face CDA may intersect with the ray
        uDC = plRay * Pluecker(vert[3], vert[2]);
        signDC = glm::sign(uDC);

        int signCDA = signAC;
        if (signCDA == 0)  {
            signCDA = -signAD;
            if (signCDA == 0)  {
                signCDA = -signDC;
            }
        }

        if ((signCDA != 0)  &&
            ((-signDC == signCDA)  ||  (signDC == 0)))  {
            // Face CDA intersects with the ray
            // Face DCB also intersects with the ray
            double invVolCDA = 1.0 / (uAC - uDC - uAD);

            if (signBC == -2)  {
                uBC = plRay * Pluecker(vert[1], vert[2]);
            }
            if (signDB == -2)  {
                uDB = plRay * Pluecker(vert[3], vert[1]);
            }
            double invVolDCB = 1.0 / (uDC - uBC - uDB);
            
            if (signCDA == 1)  {
                enterFace = 1;
                uEnter1 =  uAC * invVolCDA;
                uEnter2 = -uDC * invVolCDA;
                enterPoint =   (1-uEnter1-uEnter2)*vert[2]
                             + uEnter1*vert[3] + uEnter2*vert[0];
                
                leaveFace = 0;
                uLeave1 = -uDB * invVolDCB;
                uLeave2 =  uDC * invVolDCB;
                leavePoint =   (1-uLeave1-uLeave2)*vert[3]
                             + uLeave1*vert[2] + uLeave2*vert[1];

                ComputeParametricDist(orig, dir,
                                      enterPoint, leavePoint,
                                      tEnter, tLeave);
                return true;
            }
            else  {
                leaveFace = 1;
                uLeave1 =  uAC * invVolCDA;
                uLeave2 = -uDC * invVolCDA;
                leavePoint =   (1-uLeave1-uLeave2)*vert[2]
                             + uLeave1*vert[3] + uLeave2*vert[0];
                
                enterFace = 0;
                uEnter1 = -uDB * invVolDCB;
                uEnter2 =  uDC * invVolDCB;
                enterPoint =   (1-uEnter1-uEnter2)*vert[3] 
                             + uEnter1*vert[2] + uEnter2*vert[1];

                ComputeParametricDist(orig, dir,
                                      enterPoint, leavePoint,
                                      tEnter, tLeave);
                return true;
            }
        }
    }

    // Three faces do not intersect with the ray, the fourth will not.
    return false;
}



