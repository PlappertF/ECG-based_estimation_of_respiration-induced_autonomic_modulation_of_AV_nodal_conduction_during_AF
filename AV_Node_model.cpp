#include "mex.h"
#include <math.h>
#include <vector>
#include <queue>

using namespace std;

inline double compR(double DI,double *p, double t, double *a)  { return (1+a[0] * sin( a[1] * 2 * M_PI * t / 1000)/2)*(p[0]+p[1]*(1.0-exp(-DI/p[2]))); }
inline double compD(double DI,double *p, double t, double *a)  { return (1+a[0] * sin( a[1] * 2 * M_PI * t / 1000)/2)*(p[0]+p[1]*exp(-DI/p[2])); }

struct imp
{
    double time; // Time at which the impulse will arrive at the node
    int node; // Node at which the impulse will arrive
    int imp_idx; // index of impulse which enters the AV node
    int path; // path on which the impulse started, this is interesting for
              // AVN reentry, when the excitation from one pathway returns 
              // on the other pathway
};

class cmpimp {
public:
    bool operator()(imp& i1, imp& i2)
    {
        if (i1.time > i2.time) return true;
        return false;
    }
};

// The arguments nlhs and nrhs contain the number of left side and right 
// side arguments, respectively. prhs is an array of mxArray pointers whose
// length is nrhs. plhs is an array whose length is nlhs, where your
// function must set pointers for the output mxArrays.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // nlhs - Number of output arguments (int): Number of expected mxArray output
    //        arguments, specified as an integer.
    // plhs - MATLAB arrays (mxArray*): Array of pointers to the expected
    //        mxArray output arguments.
    // nrhs - Number of input arguments (int): Number of input mxArrays,
    //        specified as an integer.
    // prhs - MATLAB arrays (const mxArray*): Array of pointers to the 
    //        mxArray input arguments. Do not modify any prhs values in 
    //        your MEX file. Changing the data in these read-only mxArrays
    //        can produce undesired side effects.
    
// INPUT    
    double *AVEntries =   (double *) mxGetPr(prhs[0]); // Series of atrial arrival times
    double *apd_f =     (double *) mxGetPr(prhs[1]); // Refractory period parameters of fast pathway
    double *apd_s =     (double *) mxGetPr(prhs[2]); // Refractory period parameters of slow pathway
    double *delay_f =   (double *) mxGetPr(prhs[3]); // Conduction delay parameters of fast pathway
    double *delay_s =   (double *) mxGetPr(prhs[4]); // Conduction delay parameters of slow pathway
    double *A_R =       (double *) mxGetPr(prhs[5]); // ANS factor that affects the refractory period

    int nAVEntries = mxGetNumberOfElements(prhs[0]);

// OUTPUT
    plhs[0] = mxCreateDoubleMatrix(nAVEntries,2,mxREAL); // Series of ventricular arrival times
    // HisEntries has two columns to track the ventricular activations for the slow and fast pathway separately for the cases that the same atrial impulse is causing two ventricular impulses
    double* HisEntries = mxGetPr(plhs[0]);
    for (int i=0; i<2*nAVEntries; i++) HisEntries[i] = mxGetNaN();
    
    plhs[1] = mxCreateDoubleMatrix(nAVEntries,1,mxREAL); 
    // For each node, it is stored if the node was excited by a wave 
    // originating in the slow or fast pathway. In the case of a concealed
    // conduction it can be analysed if it is because of a long refractory
    // period from the previous excitation, or because of a returning 
    // wave from AV node reentry.
    // For each atrial impulse, this array stores whether this impulse resulted in a ventricular excitation and over which pathway the impulse traveled.
    // 0: Concealed conduction
    // 1: The atrial impulse resulted in a ventricular excitation and travelled over the slow pathway
    // 2: The atrial impulse resulted in a ventricular excitation and travelled over the fast pathway
    // 3: The atrial impulse resulted in two ventricular excitation, where one excitation resulted from the impulse travelling over the slow pathway and the other from an impulse from the slow pathway
    double* AHPathway = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(nAVEntries,22,mxREAL); // This matrix tracks the refractory periods of all nodes and all atrial impulses
    double* AHRT = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(nAVEntries,22,mxREAL); // This matrix tracks the conduction delay of all nodes and all atrial impulses
    double* AHDL = mxGetPr(plhs[3]);
    
    plhs[4] = mxCreateDoubleMatrix(nAVEntries,22,mxREAL); // This matrix tracks at which time a node was excited, listed for all nodes and all atrial impulses
    double* AHET = mxGetPr(plhs[4]);
    
    plhs[5] = mxCreateDoubleMatrix(22,1,mxREAL);
    // last refractory time value of all 22 nodes, The coupling node is modelled as two separate node for slow and fast pathway to track that one incoming excitation can cause two outgoing excitations
    double* RT = mxGetPr(plhs[5]);
    // This for loop is only one line long and fills the RT with 21 zeros.
    for (int i=0; i<22; i++) RT[i] = 0;
    
    plhs[6] = mxCreateDoubleMatrix(22,1,mxREAL);
    // last conduction delay value of all 22 nodes,
    double* DL = mxGetPr(plhs[6]);
    
    plhs[7] = mxCreateDoubleMatrix(22,1,mxREAL);
    // Tracks for each node, if the impulse originated in the slow or fast pathway
    double* PA = mxGetPr(plhs[7]);
    
    //int* index = (int *) mxCalloc(21, sizeof(int));
    int current_index(0);
    
    int i(0), counter(0);
    double *next_apd, *next_delay, *next_AR;

    double end_refper;
    int next_node;
    
    // Create a queue q. For every incoming excitation wave, two entries
    // are created, once for the fast and slow pathway, respectively. The 
    // struct s is filled for every excitation wave entering both pathways.
    // When the struct s is full, it is pushed in the queue and the
    // struct s is filled for the next entry.
    priority_queue<imp, vector<imp>, cmpimp> q;
    imp s;
    for (int i = 0; i < nAVEntries; ++i)
    {   
        // time and imp_idx is initialized for both the fast and slow path
        s.time = AVEntries[i];
        s.imp_idx = i;
        // slow pathway going from node 0 to node 9
        s.node = 0;
        s.path = 1;
        q.push(s);
        // fast pathway going from node 10 to node 19
        //s.time = AVEntries[i]+20;
        s.node = 10;
        s.path = 2;
        q.push(s);
    }
    
    // This while loop is running until all excitation waves are processed
    // in the queue q.   
    while (!q.empty())
    {     
        s = q.top();
        q.pop();
                
        next_node = s.node;
        if (next_node<10) // slow pathway
        {
            next_apd = apd_s;
            next_delay = delay_s;
        }
        else if (next_node<20) // fast pathway
        {
            next_apd = apd_f;
            next_delay = delay_f;
        }
        next_AR = A_R;
        
        if ( s.time >= RT[next_node] ) // Loop is executed if the refractory period of the current node is not ongoing
        {         
            if (next_node<20) // Only calculate this for the slow and fast pathway
            {
                // Computes conduction delay for the current impulse to neighboring nodes
                DL[next_node] = compD(s.time-RT[next_node],next_delay,s.time,next_AR);
                // Computes refractory period for the current node due to the current excitation
                RT[next_node] = s.time + compR(s.time-RT[next_node],next_apd,s.time,next_AR);
            }
            else // For End node
            {
                DL[next_node] = 0;
                RT[next_node] = s.time + 250;  
            }

            // Stores the refractory period, conduction delay and excitation time in the matrix to analyse the model behavior afterwards in Matlab
            AHRT[s.imp_idx+next_node*nAVEntries] = RT[next_node];
            AHDL[s.imp_idx+next_node*nAVEntries] = DL[next_node];
            AHET[s.imp_idx+next_node*nAVEntries] = s.time;

            if (next_node == 0 || next_node == 10) // first node of slow OR fast pathway
            {
                s.time = s.time+DL[next_node];
                s.node = next_node+1; // impulse moves forwards
                q.push(s);
            }
            else if (next_node == 9) // last node of slow pathway
            {
                s.time = s.time+DL[next_node];
                s.node = 8; // impulse moves backwards on slow pathway
                q.push(s);
                // The single coupling node is modelled as two nodes. 
                // The two nodes are synchronized to function like a single node. 
                // Modelling the node as two allows for tracking of the 
                // refractory period, conduction delay and excitation time 
                // for impulses that resulted in two ventricular excitations 
                // and therefore passed the coupling node twice.
                if (s.path == 1) // impulse originated in slow pathway
                {
                    s.node = 20; // impulse moves to end node for slow pathway impulses
                }
                else // s.path == 2  impulse originated in fast pathway
                {
                    s.node = 21; // impulse moves to end node for fast pathway impulses
                }
                q.push(s);
                s.node = 19; //impulse moves to last node of fast pathway
                q.push(s);                  
            }
            else if (next_node == 19)
            {
                s.time = s.time+DL[next_node];
                s.node = 18; // impulse moves backwards on fast pathway
                q.push(s);
                if (s.path == 1) // impulse originated in slow pathway
                {
                    s.node = 20; // impulse moves to end node for slow pathway impulses
                }
                else // s.path == 2  impulse originated in fast pathway
                {
                    s.node = 21; // impulse moves to end node for fast pathway impulses
                }
                q.push(s);
                s.node = 9; // impulse moves to last node of slow pathway
                q.push(s);                     
            }  
            else if (next_node<19) // center nodes of slow OR fast pathway
            {
                s.time = s.time+DL[next_node];
                s.node = next_node-1; // impulse moves backwards
                q.push(s);
                s.node = next_node+1; // impulse moves forwards
                q.push(s);
            }
            else if (next_node == 20) // end node for slow pathway impulses
            {                
                HisEntries[s.imp_idx] = s.time+DL[next_node];
                AHPathway[s.imp_idx] = AHPathway[s.imp_idx] + s.path;
                DL[21] = DL[next_node]; // synchronize the conduction delay of the two end node, since they are technically one single node
                RT[21] = RT[next_node]; // synchronize the refractory period of the two end node, since they are technically one single node
            }
            else if (next_node == 21) // end node for fast pathway impulses
            {                
                HisEntries[s.imp_idx+nAVEntries] = s.time+DL[next_node];
                AHPathway[s.imp_idx] = AHPathway[s.imp_idx] + s.path;
                DL[20] = DL[next_node]; // synchronize the conduction delay of the two end node, since they are technically one single node
                RT[20] = RT[next_node]; // synchronize the refractory period of the two end node, since they are technically one single node
            }
        
        } // end if  
        
    } // end while
    //mxFree(DL);

    //mxFree(index);
}
