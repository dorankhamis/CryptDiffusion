import numpy as np

def Poly_interp4(y, x, xx, ii):
   ## args: double complex[::1] y, double x, double[::1] xx, int ii
   ''' Interpolate the value of a vector y (defined on a grid xx) 
       at a point x, passing through the points xx[ii] -- x[ii+5]. '''
   a0 = (x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii]-xx[ii+1])*(xx[ii]-xx[ii+2])*(xx[ii]-xx[ii+3])*(xx[ii]-xx[ii+4])) *y[ii]
   a1 = (x-xx[ii])*(x-xx[ii+2])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii+1]-xx[ii])*(xx[ii+1]-xx[ii+2])*(xx[ii+1]-xx[ii+3])*(xx[ii+1]-xx[ii+4])) *y[ii+1]
   a2 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii+2]-xx[ii])*(xx[ii+2]-xx[ii+1])*(xx[ii+2]-xx[ii+3])*(xx[ii+2]-xx[ii+4])) *y[ii+2]
   a3 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+4]) / ((xx[ii+3]-xx[ii])*(xx[ii+3]-xx[ii+1])*(xx[ii+3]-xx[ii+2])*(xx[ii+3]-xx[ii+4])) *y[ii+3]
   a4 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+3]) / ((xx[ii+4]-xx[ii])*(xx[ii+4]-xx[ii+1])*(xx[ii+4]-xx[ii+2])*(xx[ii+4]-xx[ii+3])) *y[ii+4]
   return a0+a1+a2+a3+a4

def integrate(f, x):
   ## args: double complex[::1] f, double[::1] x
   ''' Integrate an array f on a grid x '''
   x = np.asarray(x)
   f = np.asarray(f)
   n = x.shape[0]
   s = np.zeros(n, dtype=np.float32)
   for i in range(1,n-5):
      jj=i-1 
      bma = x[i]-x[i-1]
      x1 = x[i-1] + bma/4.; x2 = x[i-1] + bma/2.; x3 = x[i-1] + 3.*bma/4.
      s[i] = s[i-1] + bma/90. *(7.*f[i-1] + 32.*Poly_interp4(f, x1, x, jj) + 12.*Poly_interp4(f, x2, x, jj) + 32.*Poly_interp4(f, x3, x, jj) + 7.*f[i])
   jj = n-5
   for i in range(n-5,n):
      bma = x[i]-x[i-1]
      x1 = x[i-1] + bma/4.; x2 = x[i-1] + bma/2.; x3 = x[i-1] + 3.*bma/4.
      s[i] = s[i-1] + bma/90. *(7.*f[i-1] + 32.*Poly_interp4(f, x1, x, jj) + 12.*Poly_interp4(f, x2, x, jj) + 32.*Poly_interp4(f, x3, x, jj) + 7.*f[i])
   return s
