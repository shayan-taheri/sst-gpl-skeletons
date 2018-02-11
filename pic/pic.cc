#include <vector>
#include <cstdlib>
#include <mpi.h>

struct Particle {
  double x[3];
  double v[3];
  double m;
  double deltaT;
  int cell;
};

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
#pragma sst null_type sstmac::vector size resize 
  std::vector<Particle> local; 
  std::vector<Cell> cells;
  std::vector<Migration> outgoing;
  std::vector<Migration> incoming;
  int od[3]; //the od factor in each dim
  int localGridDims[3];
  double microIterScale[4][4][4];
  double migrateFraction[4][4][4];
  int boxOcc[4][4][4]; //allow for max 4x overdecomp in each dim
};

void skeletonInitOutgoing(Patch& p){
  for (int x=0; x < p.od[0]; x++){
    for (int y=0; y < p.od[1]; y++){
      for (int z=0; z < p.od[2]; z++){
        p.migrateFraction[x][y][z] = p.microIterScale[x][y][z];
      }
    }
  }
}

int skeletonFillOutgoing(Patch& p){
  int nbr = 0;
  int numMoving = 0;
  int totalMoving = 0;
  for (int y=0; y < p.od[0]; ++y){
    for (int z=0; z < p.od[2]; ++z){
      //these go left on x=0
      double frac = p.migrateFraction[0][y][z];
      numMoving += frac * p.boxOcc[0][y][z];
      p.migrateFraction[0][y][z] *= p.microIterScale[0][y][z];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  ++nbr;
  numMoving = 0;
  int lastX = p.od[0] - 1;
  for (int y=0; y < p.od[1]; ++y){
    for (int z=0; z < p.localGridDims[2]; ++z){
      //these go left on x=0
      double frac = p.migrateFraction[lastX][y][z];
      numMoving += frac * p.boxOcc[lastX][y][z];
      p.migrateFraction[lastX][y][z] *= p.microIterScale[lastX][y][z];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  ++nbr;
  numMoving = 0;
  for (int x=0; x < p.od[0]; ++x){
    for (int z=0; z < p.od[2]; ++z){
      //these go left on x=0
      double frac = p.migrateFraction[x][0][z];
      numMoving += frac * p.boxOcc[x][0][z];
      p.migrateFraction[x][0][z] *= p.microIterScale[x][0][z];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  ++nbr;
  numMoving = 0;
  int lastY = p.od[1] - 1;
  for (int x=0; x < p.od[0]; ++x){
    for (int z=0; z < p.od[2]; ++z){
      //these go left on x=0
      double frac = p.migrateFraction[x][lastY][z];
      numMoving += frac * p.boxOcc[x][lastY][z];
      p.migrateFraction[x][lastY][z] *= p.microIterScale[x][lastY][z];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  ++nbr;
  numMoving = 0;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.od[1]; ++y){
      //these go left on x=0
      double frac = p.migrateFraction[x][y][0];
      numMoving += frac * p.boxOcc[x][y][0];
      p.migrateFraction[x][y][0] *= p.microIterScale[x][y][0];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  ++nbr;
  numMoving = 0;
  int lastZ = p.od[2] - 1;
  for (int x=0; x < p.od[0]; ++x){
    for (int y=0; y < p.localGridDims[1]; ++y){
      //these go left on x=0
      double frac = p.migrateFraction[x][y][lastZ];
      numMoving += frac * p.boxOcc[x][y][lastZ];
      p.migrateFraction[x][y][lastZ] *= p.microIterScale[x][y][lastZ];
    }
  }
  p.outgoing[nbr].parts.resize(numMoving);
  totalMoving += numMoving;

  return totalMoving;
}

void moveParticle(Particle& part, Patch& patch)
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
    }
  }
}

int exchange(Patch& patch)
{
  std::vector<int> numSending(patch.outgoing.size());
  std::vector<int> numRecving(patch.incoming.size());
  std::vector<MPI_Request> sizeRequests(patch.outgoing.size() + patch.incoming.size());
  std::vector<MPI_Request> partRequests;
  MPI_Request collectiveReq;
  int totalOutgoing = 42; //to get the while loop going
  int systemTotal = 0;
  while (systemTotal > 0){
    totalOutgoing = 0;
    for (int f=0; f < patch.outgoing.size(); ++f){
      Migration& m = patch.outgoing[f];
      numSending[f] = m.parts.size();
      totalOutgoing += m.parts.size();
      MPI_Isend(&numSending[f], 1, MPI_INT, m.rank, send_size_tag,
                MPI_COMM_WORLD, &sizeRequests[f]);
      if (numSending[f] > 0){
        MPI_Request req;
        MPI_Isend(m.parts.data(), m.parts.size() * sizeof(Particle), MPI_BYTE,
                  m.rank, send_parts_tag, MPI_COMM_WORLD, &req);
        partRequests.push_back(req);
      }
    }

    for (int f=0; f < patch.incoming.size(); ++f){
      Migration& m = patch.incoming[f];
      MPI_Irecv(&numRecving[f], 1, MPI_INT, m.rank, send_size_tag,
                MPI_COMM_WORLD, &sizeRequests[f + patch.outgoing.size()]);
    }

    MPI_Iallreduce(&totalOutgoing, &systemTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &collectiveReq);

    int numDone = 0;
    int numNeeded = patch.incoming.size() + patch.outgoing.size();
    while (numDone < numNeeded){
      int idx;
      MPI_Waitany(numNeeded, sizeRequests.data(), &idx, MPI_STATUSES_IGNORE);
      if (idx > patch.outgoing.size()){
        int recvIdx = idx - patch.outgoing.size();
        //we received the size of a particle we need - post its recv
        int numIncoming = numRecving[recvIdx];
        if (numIncoming > 0){
          Migration& m = patch.incoming[recvIdx];
          m.parts.resize(numIncoming);
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
  }
  return systemTotal;
}

void move(Patch& patch)
{
#pragma omp parallel for
  for (int i=0; i < patch.local.size(); ++i){
    Particle& part = patch.local[i];
    part.deltaT = patch.deltaT;
    moveParticle(part, patch);
  }

#pragma sst call skeletonInitOutgoing(patch)
#pragma sst init skeletonFillOutgoing(patch)
  int systemTotalMoves = exchange(patch);

  while (systemTotalMoves > 0){
#pragma omp parallel for
    for (int i=0; i < patch.local.size(); ++i){
      Particle& part = patch.local[i];
      moveParticle(part, patch);
    }
#pragma sst call skeletonFillOutgoing(patch)
    exchange(patch);
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
  //id = z*ny*nx + y*nx + x;
  int myZ = patch.id / (nPatchesX*nPatchesY);
  int remId = patch.id % (nPatchesX*nPatchesY);
  int myY = remId / nPatchesX;
  int myX = remId % nPatchesX;

  patch.outgoing.resize(6);
  patch.incoming.resize(6);

  //I have a plus X partner
  int plusX = (myX + 1) % nPatchesX;
  int plusXpartner = patchId(plusX,myY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[0].rank = plusXpartner;
  patch.incoming[1].rank = plusXpartner;


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

  //I have a minus X partner
  int minusX = (myX + nPatchesX - 1) % nPatchesX;
  int minusXpartner = patchId(minusX,myY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[1].rank = minusXpartner;
  patch.incoming[0].rank = minusXpartner;

  //I have a plus Y partner
  int plusY = (myY + 1) % nPatchesY;
  int plusYpartner = patchId(myX,plusY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[2].rank = plusYpartner;
  patch.incoming[3].rank = plusYpartner;

  //I have a minus Y partner
  int minusY = (myY + nPatchesY - 1) % nPatchesY;
  int minusYpartner = patchId(myX,minusY,myZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[3].rank = minusYpartner;
  patch.incoming[2].rank = minusYpartner;

  //I have a plus Z partner
  int plusZ = (myZ + 1) % nPatchesZ;
  int plusZpartner = patchId(myX,myY,plusZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[4].rank = plusZpartner;
  patch.incoming[5].rank = plusZpartner;

  //I have a minus Z partner
  int minusZ = (myZ + nPatchesZ - 1) % nPatchesZ;
  int minusZpartner = patchId(myX,myY,minusZ,nPatchesX,nPatchesY,nPatchesZ);
  patch.outgoing[5].rank = minusZpartner;
  patch.incoming[4].rank = minusZpartner;
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
  myPatch.od[0] = 2;
  myPatch.od[1] = 2;
  myPatch.od[2] = 2;
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

  init(myPatch, ppc, nPatchesX, nPatchesY, nPatchesZ);
  for (int s=0; s < nSteps; ++s){
    move(myPatch);
  }

  MPI_Finalize();
  return 0;
}

