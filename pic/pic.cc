#include <vector>
#include <cstdlib>
#include <mpi.h>
#include <sprockit/keyword_registration.h>

RegisterKeywords(
 {"odx", "the level of intra-patch overdecomposition in X direction"},
 {"ody", "the level of intra-patch overdecomposition in Y direction"},
 {"odz", "the level of intra-patch overdecomposition in Z direction"},
 {"unitx", "the stride of X units with uniform particle move characteristics" },
 {"unity", "the stride of Y units with uniform particle move characteristics" },
 {"unitz", "the stride of Z units with uniform particle move characteristics" },
 {"scramble", "whether to constantly scramble particle movements or keep steady flow"},
 {"max_migration", "the max fraction of particles that can migrate"},
 {"min_migration", "the min fraction of particles that can migrate"}
);

struct Particle {
  double x[3];
  double v[3];
  double m;
  double deltaT;
  int cell;
};

//0th incoming face is +X
char incomingChars[] = { '+', '-' };
//0th outgoing face is -X
char outgoingChars[] = { '-', '+' };
char faceChars[] = {'X', 'Y', 'Z' };

#define inChar(x) incomingChars[x%2]
#define outChar(x) outgoingChars[x%2]
#define dimChar(x) faceChars[x/2]

#define debug(...) //printf(__VA_ARGS__)

#define SKEL_MAX_OD 4
#define SKEL_NUM_FACES 6

static const int send_size_tag = 100;
static const int send_parts_tag = 101;

struct Migration {
#pragma sst null_type sstmac::vector size resize clear data
  std::vector<Particle> parts;
  int rank; 
};

double dot(double a[3], double b[3]){
  double prod = 0;
  prod += a[0]*b[0];
  prod += a[1]*b[1];
  prod += a[2]*b[2];
  return prod;
}

struct Face {
  double n[3];
  double x[3];
  int dstCell;
  int dstRank;
};

struct Cell {
  std::vector<Face> faces;
};

struct Patch {
  int id;
  int nPatches;
  int nCells;
  double center[3];
  double spacing;
  double deltaT;
  int localGridDims[3];
  int gridPosition[3];
#pragma sst null_type sstmac::vector size resize 
  std::vector<Particle> local; 
#pragma sst null_type sstmac::vector size resize push_back
  std::vector<int> holes;
#pragma sst null_type sstmac::vector size resize
  std::vector<Cell> cells;
  std::vector<Migration> outgoing;
  std::vector<Migration> incoming;

  /** The extra variables for running the skeleton */
  int od[3]; //the od factor in each dim
  int uniformityFactory[3];
  //What decrease fraction per micro-iteration in # particles moving
  double microIterScale;
  //What is current fraction of particles migrating
  double migrateFraction;
  double minMigrateFraction;
  double maxMigrateDifference;
  bool scrambleMigration;
  int boxOcc[SKEL_MAX_OD][SKEL_MAX_OD][SKEL_MAX_OD]; //allow for max 4x overdecomp in each dim
};

static inline double getSkeletonFraction(int step, Patch& p, int dim)
{
  static const int primes[] = {
   1571, 2459, 2801, 3559, 3079, 3019, 6269
  };
  static const int numPrimes = sizeof(primes) / sizeof(int);

  int myIdx = p.od[dim]*p.gridPosition[dim] / p.uniformityFactory[dim];
  if (p.scrambleMigration){
    myIdx = primes[(myIdx+step) % numPrimes] % numPrimes;
  }
  double myExtraFraction = p.maxMigrateDifference * myIdx / numPrimes;
  return p.minMigrateFraction + myExtraFraction;
}

void skeletonInitOutgoing(int step, Patch& p)
{
  double xInc = getSkeletonFraction(step, p, 0);
  double yInc = getSkeletonFraction(step, p, 1);
  double zInc = getSkeletonFraction(step, p, 2);

  p.microIterScale = p.migrateFraction = (xInc + yInc + zInc);

  debug("Rank %d set scale factor to %f\n", p.id, p.microIterScale);
}



void skeletonInitOverdecomposition(Patch& p, int ppc){

  auto params = get_params();

  p.od[0] = params->get_optional_int_param("odx", 1);
  p.od[1] = params->get_optional_int_param("ody", 1);
  p.od[2] = params->get_optional_int_param("odz", 1);
  p.uniformityFactory[0] = params->get_optional_int_param("unitx", 1);
  p.uniformityFactory[1] = params->get_optional_int_param("unity", 1);
  p.uniformityFactory[2] = params->get_optional_int_param("unitz", 1);
  int numOdBoxes = p.od[0] * p.od[1] * p.od[2];
  int numPartsPerOdBox = (ppc*p.nCells) / numOdBoxes;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        p.boxOcc[x][y][z] = numPartsPerOdBox;
      }
    }
  }

  p.minMigrateFraction = params->get_optional_double_param("min_migration", 0.02);
  double maxMig = params->get_optional_double_param("max_migration", 0.06);
  p.maxMigrateDifference = maxMig - p.minMigrateFraction;
  if (p.maxMigrateDifference < 0){
    std::cerr << "Max migration=" << maxMig << " must be greater than min="
              << p.minMigrateFraction << std::endl;
    abort();
  }
  p.scrambleMigration = params->get_optional_bool_param("scramble", false);
}

void skeletonPackMigrated(Patch& p){
  int totalIncoming = 0;
  { Migration& m = p.incoming[0];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[1] * p.od[2]);
    totalIncoming += m.parts.size();
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        p.boxOcc[0][y][z] += numIncomingPerBox;
        debug("Rank %d box %d-%d-%d now at %d parts after its share of %d\n",
               p.id,0,y,z,p.boxOcc[0][y][z],int(m.parts.size()));
      }
    }
    m.parts.clear();
  } }

  { Migration& m = p.incoming[1];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[1] * p.od[2]);
    totalIncoming += m.parts.size();
    int lastX = p.od[0] - 1;
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        p.boxOcc[lastX][y][z] += numIncomingPerBox;
      }
    }
    m.parts.clear();
  } }

  { Migration& m = p.incoming[2];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[0] * p.od[2]);
    totalIncoming += m.parts.size();
    for (int x=0; x < p.od[0]; ++x){
      for (int z=0; z < p.od[2]; ++z){
        p.boxOcc[x][0][z] += numIncomingPerBox;
      }
    }
    m.parts.clear();
  } }

  { Migration& m = p.incoming[3];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[0] * p.od[2]);
    totalIncoming += m.parts.size();
    int lastY = p.od[1] - 1;
    for (int x=0; x < p.od[0]; ++x){
      for (int z=0; z < p.od[2]; ++z){
        p.boxOcc[x][lastY][z] += numIncomingPerBox;
      }
    }
    m.parts.clear();
  } }

  { Migration& m = p.incoming[4];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[0] * p.od[1]);
    totalIncoming += m.parts.size();
    for (int x=0; x < p.od[0]; ++x){
      for (int y=0; y < p.od[1]; ++y){
        p.boxOcc[x][y][0] += numIncomingPerBox;
      }
    }
    m.parts.clear();
  } }

  { Migration& m = p.incoming[5];
  if (m.parts.size() > 0){
    //spread evenly across the boxes
    int numIncomingPerBox = m.parts.size() / (p.od[0] * p.od[1]);
    totalIncoming += m.parts.size();
    int lastZ = p.od[2] - 1;
    for (int x=0; x < p.od[0]; ++x){
      for (int y=0; y < p.od[1]; ++y){
        p.boxOcc[x][y][lastZ] += numIncomingPerBox;
      }
    }
    m.parts.clear();
  } }

  size_t oldSize = p.local.size();
  size_t newSize = oldSize + totalIncoming;
  p.local.resize(newSize);

  int boxSum = 0;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        boxSum += p.boxOcc[x][y][z];
        debug("Rank %d box %d-%d-%d now at %d parts\n",
               p.id,x,y,z,p.boxOcc[x][y][z]);
      }
    }
  }
  debug("Rank %d sum over boxes is %d\n", p.id, boxSum);
}


int skeletonFillOutgoing(Patch& p){
  int occIncrease[SKEL_MAX_OD][SKEL_MAX_OD][SKEL_MAX_OD];
  int numMovingPerFace[SKEL_NUM_FACES];

  ::memset(occIncrease, 0, sizeof(occIncrease));
  ::memset(numMovingPerFace, 0, sizeof(numMovingPerFace));

  int lastX = p.od[0] - 1;
  int lastY = p.od[1] - 1;
  int lastZ = p.od[2] - 1;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        double frac = p.migrateFraction;
        int boxMoves = frac * p.boxOcc[x][y][z];
        if (x==0 && x == lastX){
          numMovingPerFace[0] += boxMoves;
          numMovingPerFace[1] += boxMoves;
        } else if (x == 0) {
          numMovingPerFace[0] += boxMoves;
          occIncrease[x+1][y][z] += boxMoves;
        } else if (x == lastX) {
          numMovingPerFace[1] += boxMoves;
          occIncrease[x-1][y][z] += boxMoves;
        } else {
          occIncrease[x+1][y][z] += boxMoves;
          occIncrease[x-1][y][z] += boxMoves;
        }

        if (y==0 && y == lastY){
          numMovingPerFace[2] += boxMoves;
          numMovingPerFace[3] += boxMoves;
        } else if (y == 0) {
          numMovingPerFace[2] += boxMoves;
          occIncrease[x][y+1][z] += boxMoves;
        } else if (y == lastY) {
          numMovingPerFace[3] += boxMoves;
          occIncrease[x][y-1][z] += boxMoves;
        } else {
          occIncrease[x][y+1][z] += boxMoves;
          occIncrease[x][y-1][z] += boxMoves;
        }

        if (z==0 && z == lastZ){
          numMovingPerFace[4] += boxMoves;
          numMovingPerFace[5] += boxMoves;
        } else if (z == 0) {
          numMovingPerFace[4] += boxMoves;
          occIncrease[x][y][z+1] += boxMoves;
        } else if (z == lastZ) {
          numMovingPerFace[5] += boxMoves;
          occIncrease[x][y][z-1] += boxMoves;
        } else {
          occIncrease[x][y][z+1] += boxMoves;
          occIncrease[x][y][z-1] += boxMoves;
        }
        //we lose some out each face
        p.boxOcc[x][y][z] -= 6*boxMoves;
      }
    }
  }

  p.migrateFraction *= p.microIterScale;

  int totalLocalOcc = 0;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.od[1]; ++y){
      for (int z=0; z < p.od[2]; ++z){
        int oldOcc = p.boxOcc[x][y][z];
        p.boxOcc[x][y][z] += occIncrease[x][y][z];
        debug("Rank %d box %d-%d-%d depleted to %d, but back up to %d\n",
               p.id, x,y,z, oldOcc, p.boxOcc[x][y][z]);
        totalLocalOcc += p.boxOcc[x][y][z];
      }
    }
  }

  int totalMoving = 0;
  for (int i=0; i < SKEL_NUM_FACES; ++i){
    debug("Rank %d produced %d outgoing on face %c%c\n",
           p.id, numMovingPerFace[i], outChar(i), dimChar(i));
    p.outgoing[i].parts.resize(numMovingPerFace[i]);
    totalMoving += numMovingPerFace[i];
  }

  debug("Total is %d + %d = %d -> %d per\n",
         totalLocalOcc, totalMoving, totalLocalOcc + totalMoving,
         (totalLocalOcc + totalMoving) / (p.od[0]*p.od[1]*p.od[2]));

  int newSize = p.local.size() - totalMoving;
  if (totalMoving > int(p.local.size())){
    std::cerr << "more particles=" << totalMoving
              << " moving than exist=" << int(p.local.size())
              << std::endl;
    abort();
  }
  p.local.resize(newSize);

  return totalMoving;
}

void backfill(Patch& patch)
{
  //stub - to be filled in
}

//this is never actually usable
#pragma sst delete
void moveParticle(int idx, Particle& part, Patch& patch)
{
#pragma sst loop_count 2
  while (part.deltaT > 0){
    int minFace = 0;
    double minDeltaT = part.deltaT;
    //find intersection with each face
    Cell& cell = patch.cells[part.cell];
#pragma sst loop_count 6 //6 faces
    for (int f=0; f < cell.faces.size(); ++f){
      Face& face = cell.faces[f];
      double vecComponent = dot(face.n, part.v);
#pragma sst branch_predict 0.5 //this will be true on half of the faces
      if (vecComponent > 0){
        //great - will hit this face
        double deltaX[3];
        deltaX[0] = face.x[0] - part.x[0];
        deltaX[1] = face.x[1] - part.x[1];
        deltaX[2] = face.x[2] - part.x[2];
        double distance = dot(deltaX, face.n);
        double tIntersect = distance / vecComponent;
        double tMove = std::min(tIntersect, part.deltaT);
        if (tMove < minDeltaT){
          minFace = f;
          minDeltaT = tMove;
        }
      }
    }
    Face& dstFace = cell.faces[minFace];
    part.cell = dstFace.dstCell;
    if (dstFace.dstRank >= 0){
      patch.outgoing[dstFace.dstRank].parts.push_back(part);
      patch.holes.push_back(idx);
    }
  }
}

void packMigrated(Patch& patch);

int exchange(Patch& patch)
{
  std::vector<int> numSending(patch.outgoing.size(), 0);
  std::vector<int> numRecving(patch.incoming.size(), 0);
  std::vector<MPI_Request> sizeRequests(patch.outgoing.size() + patch.incoming.size());
  std::vector<MPI_Request> partRequests;
  MPI_Request collectiveReq;
  int totalOutgoing = 0;
  int systemTotal = 0;
  for (int f=0; f < patch.outgoing.size(); ++f){
    Migration& m = patch.outgoing[f];
    numSending[f] = m.parts.size();
    totalOutgoing += m.parts.size();
#pragma sst keep
    MPI_Isend(&numSending[f], 1, MPI_INT, m.rank, send_size_tag,
              MPI_COMM_WORLD, &sizeRequests[f]);
    debug("Rank %d sending %d to %d on face %c%c\n",
           patch.id, numSending[f], m.rank, outChar(f), dimChar(f));
    if (numSending[f] > 0){
      MPI_Request req;
      MPI_Isend(m.parts.data(), m.parts.size() * sizeof(Particle), MPI_BYTE,
                m.rank, send_parts_tag, MPI_COMM_WORLD, &req);
      partRequests.push_back(req);
    }
  }

  for (int f=0; f < patch.incoming.size(); ++f){
    Migration& m = patch.incoming[f];
#pragma sst keep
    MPI_Irecv(&numRecving[f], 1, MPI_INT, m.rank, send_size_tag,
              MPI_COMM_WORLD, &sizeRequests[f + patch.outgoing.size()]);
  }

  MPI_Iallreduce(&totalOutgoing, &systemTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &collectiveReq);

  int numDone = 0;
  int numNeeded = patch.incoming.size() + patch.outgoing.size();
  while (numDone < numNeeded){
    int idx;
    MPI_Waitany(numNeeded, sizeRequests.data(), &idx, MPI_STATUSES_IGNORE);
    if (idx >= patch.outgoing.size()){
      int recvIdx = idx - patch.outgoing.size();
      //we received the size of a particle we need - post its recv
      int numIncoming = numRecving[recvIdx];
      if (numIncoming > 0){
        Migration& m = patch.incoming[recvIdx];
        m.parts.resize(numIncoming);
        debug("Rank %d receiving %d from %d on face %c%c\n",
               patch.id, numRecving[recvIdx], m.rank, inChar(recvIdx), dimChar(recvIdx));
        MPI_Request req;
        MPI_Irecv(m.parts.data(), numIncoming*sizeof(Particle), MPI_BYTE,
                  m.rank, send_parts_tag, MPI_COMM_WORLD, &req);
        partRequests.push_back(req);
      }
    }
    ++numDone;
  }
  if (partRequests.size()){
    //we have pushed progress forward and now all sizes have been communicated
    //now wait on all the particles themselves to get shuffled
    MPI_Waitall(partRequests.size(), partRequests.data(), MPI_STATUSES_IGNORE);
  }
  partRequests.clear(); //for the next round
  MPI_Wait(&collectiveReq, MPI_STATUS_IGNORE);

  //pack the migrated particles into the main buffer
#pragma sst instead skeletonPackMigrated(patch)
  packMigrated(patch);

  return systemTotal;
}

void move(int step, Patch& patch)
{
  printf("Rank %d moving %d particles on step %d\n",
         patch.id, int(patch.local.size()), step);
#pragma omp parallel for
  for (int i=0; i < patch.local.size(); ++i){
    Particle& part = patch.local[i];
    part.deltaT = patch.deltaT;
    moveParticle(i, part, patch);
    backfill(patch);
  }

#pragma sst call skeletonFillOutgoing(patch)
#pragma sst call skeletonInitOutgoing(step,patch) //gets called first for now
  int numQuiesced = patch.local.size();
  int systemTotalMoves = exchange(patch);

  if (patch.local.size() < numQuiesced){
    std::cerr << "how is patch size " << patch.local.size()
              << " less than numQ " << numQuiesced
              << "???" << std::endl;
    abort();
  }

  while (systemTotalMoves > 0){
#pragma omp parallel for
    for (int i=numQuiesced; i < patch.local.size(); ++i){
      Particle& part = patch.local[i];
      moveParticle(i, part, patch);
    }
#pragma sst call skeletonFillOutgoing(patch)
    systemTotalMoves = exchange(patch);
    numQuiesced = patch.local.size();
  }
}

static inline int patchId(int x, int y, int z, int nx, int ny, int nz){
  return z*nx*ny + y*nx + x;
}

static inline int cellId(Patch& p, int x, int y, int z){
  return z*p.localGridDims[0]*p.localGridDims[1] +
      y*p.localGridDims[0] + x;
}

void init(Patch& patch, int ppc, int nPatchesX, int nPatchesY, int nPatchesZ)
{
  patch.local.resize(ppc*patch.nCells);
  //id = z*ny*nx + y*nx + x;
  int myZ = patch.id / (nPatchesX*nPatchesY);
  int remId = patch.id % (nPatchesX*nPatchesY);
  int myY= remId / nPatchesX;
  int myX = remId % nPatchesX;

  patch.gridPosition[0] = myX;
  patch.gridPosition[1] = myY;
  patch.gridPosition[2] = myZ;

  /**
  This is how you would initialize if it were a real app
  int lastX = patch.localGridDims[0] - 1;
  for (int y=0; y < patch.localGridDims[1]; ++y){
    for (int z=0; z < patch.localGridDims[2]; ++z){
      int localCell = cellId(patch, lastX, y, z);
      int remoteCell = cellId(patch, 0, y, z);
    }
  }
  */

  debug("Rank %d maps to %d-%d-%d\n", patch.id, myX, myY, myZ);

  patch.outgoing.resize(6);
  patch.incoming.resize(6);

  //I have a plus X partner
  int plusX = (myX + 1) % nPatchesX;
  int plusXpartner = patchId(plusX,myY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face +X is %d\n", patch.id, plusXpartner);
  patch.outgoing[1].rank = plusXpartner;
  patch.incoming[0].rank = plusXpartner;

  //I have a minus X partner
  int minusX = (myX + nPatchesX - 1) % nPatchesX;
  int minusXpartner = patchId(minusX,myY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face -X is %d\n", patch.id, minusXpartner);
  patch.outgoing[0].rank = minusXpartner;
  patch.incoming[1].rank = minusXpartner;

  //I have a plus Y partner
  int plusY = (myY + 1) % nPatchesY;
  int plusYpartner = patchId(myX,plusY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face +Y is %d\n", patch.id, plusYpartner);
  patch.outgoing[3].rank = plusYpartner;
  patch.incoming[2].rank = plusYpartner;

  //I have a minus Y partner
  int minusY = (myY + nPatchesY - 1) % nPatchesY;
  int minusYpartner = patchId(myX,minusY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face -Y is %d\n", patch.id, minusYpartner);
  patch.outgoing[2].rank = minusYpartner;
  patch.incoming[3].rank = minusYpartner;

  //I have a plus Z partner
  int plusZ = (myZ + 1) % nPatchesZ;
  int plusZpartner = patchId(myX,myY,plusZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face +Z is %d\n", patch.id, plusZpartner);
  patch.outgoing[5].rank = plusZpartner;
  patch.incoming[4].rank = plusZpartner;

  //I have a minus Z partner
  int minusZ = (myZ + nPatchesZ - 1) % nPatchesZ;
  int minusZpartner = patchId(myX,myY,minusZ,nPatchesX,nPatchesY,nPatchesZ);
  debug("Rank %d outgoing face -Z is %d\n", patch.id, minusZpartner);
  patch.outgoing[4].rank = minusZpartner;
  patch.incoming[5].rank = minusZpartner;
}

#define crash_main(rank,...) \
  if (rank == 0){ \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    usage(); \
    return 1; \
  } else { \
    return 0; \
  }

void usage()
{
  fprintf(stderr, "usage: run <nsteps> <ppc> <cells/x> <cells/y> <cells/z> "
          "<patches/x> <patches/y> <patches/z>\n");
  fflush(stderr);
}

#define sstmac_app_name pic

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  Patch myPatch;
  MPI_Comm_rank(MPI_COMM_WORLD, &myPatch.id);
  MPI_Comm_size(MPI_COMM_WORLD, &myPatch.nPatches);

  if (argc != 9){
    crash_main(myPatch.id, "bad number of arguments");
  }



  int nSteps = atoi(argv[1]);
  int ppc = atoi(argv[2]);
  //the number of cells locally
  myPatch.localGridDims[0] = atoi(argv[3]);
  myPatch.localGridDims[1]  = atoi(argv[4]);
  myPatch.localGridDims[2]  = atoi(argv[5]);
  myPatch.nCells = myPatch.localGridDims[0] * myPatch.localGridDims[1]
                    * myPatch.localGridDims[2];
  if (myPatch.nCells == 0){
    crash_main(myPatch.id, "either got zero cell dim or misformatted number")
  }

  //the number of patches locally
  int nPatchesX = atoi(argv[6]);
  int nPatchesY = atoi(argv[7]);
  int nPatchesZ = atoi(argv[8]);

  int nPatchesTotal = nPatchesX * nPatchesY * nPatchesZ;



  if (myPatch.nPatches != nPatchesTotal){
    crash_main(myPatch.id, "requested %d=%dx%dx%d patches, but have %d MPI ranks",
               nPatchesTotal, nPatchesX, nPatchesY, nPatchesZ, myPatch.nPatches);
  }

#pragma sst call skeletonInitOverdecomposition(myPatch,ppc)
  init(myPatch, ppc, nPatchesX, nPatchesY, nPatchesZ);
  for (int s=0; s < nSteps; ++s){
    move(s,myPatch);
  }

  MPI_Finalize();
  return 0;
}

