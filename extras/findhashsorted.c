#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[])
{
    int mid, m, high, low, i;
    double *position,*hashtable,*hashsearch,hashval;
    
    hashtable  = mxGetPr(prhs[0]);
    hashsearch = mxGetPr(prhs[1]);
    hashval    = hashsearch[0];
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    position = mxGetPr(plhs[0]);
    position[0] = 0;
    
    m = mxGetNumberOfElements(prhs[0]);
    
    if (m>0)
    {
        low = 0;
        high = m-1;
        i = 0;
        if (hashtable[0] > hashval)
        {
            // No way, the number is smaller than lowest number in table
            position[0] = 0;
        }
        else if (hashtable[high] < hashval)
        {
            position[0] = 0;
        }
        else
        {
            mid = (high+low)/2;
            i = 0;
            while (mid < high & low < mid)
            {
                i++;
                if (i>=2*m)
                {
                    //mexPrintf("What!! %i",mid);
                    break;
                }
                if (hashtable[mid] < hashval)
                {
                    // mexPrintf("Low %i %i %i\n",low,mid,high);
                    low = mid;
                }
                else if (hashtable[mid] > hashval)
                {
                    // mexPrintf("High %i %i %i\n",low,mid,high);
                    high = mid;
                }
                else if  (hashtable[mid] == hashval)
                {
                    // mexPrintf("Match %i %i %i\n",low,mid,high);
                    position[0] = mid+1;
                    break;
                }
                mid = (high+low)/2;
            }
            
            if (hashtable[low] == hashval)
            {
                position[0] = low+1;
            }
            else if (hashtable[high] == hashval)
            {
                position[0] = high+1;
            }
        }
    }
    
    if (position[0] == 0)
    {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
}
