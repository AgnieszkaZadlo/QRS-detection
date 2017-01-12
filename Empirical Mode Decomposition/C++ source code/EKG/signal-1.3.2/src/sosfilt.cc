// Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France) <ecm@apc.univ-paris7.fr>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#include <octave/config.h>
#include <octave/defun-dld.h>
#include <octave/error.h>
#include <octave/gripes.h>
#include <octave/oct-obj.h>
#include <octave/pager.h>
#include <octave/quit.h>
#include <octave/variables.h>

DEFUN_DLD (sosfilt, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{y} =} sosfilt (@var{sos}, @var{x})\n\
Second order section IIR filtering of @var{x}.  The second order section\n\
filter is described by the matrix @var{sos} with:\n\
\n\
@multitable {col 1} {this is column two}\n\
@item @tab [ @var{B1} @var{A1} ]@*\n\
@item @var{sos} = @tab [ @dots{} ],@*\n\
@item @tab [ @var{BN} @var{AN} ]@*\n\
@end multitable\n\
\n\
where @code{@var{B1} = [b0 b1 b2]} and @code{@var{A1} = [1 a1 a2]} for\n\
section 1, etc.  The b0 entry must be nonzero for each section.\n\
@end deftypefn\n")
{
  octave_value_list retval;

  int nargin = args.length ();

  if (nargin != 2)
    {
      print_usage ();
      return retval;
    }

  Matrix sos( args(0).matrix_value() );

  if (error_state)
    {
      gripe_wrong_type_arg("sosfilt",args(0));
      return retval;
    }

  if (sos.columns() != 6)
    {
      error("Second-order section matrix must be a non-empty Lx6 matrix");
      return retval;
    }

  Matrix x( args(1).matrix_value() );

  if (error_state)
    {
      gripe_wrong_type_arg("sosfilt",args(1));
      return retval;
    }

  int n=x.rows();
  int m=x.columns();

  bool isrowvector=false;

  if ((n==1)&&(m>1)) // if row vector, transpose to column vector
    {
      x=x.transpose();
      n=x.rows();
      m=x.columns();
      isrowvector=true;
    }

  Matrix y(n,m,0.0);

  for (int k=0; k<m; k++)
    {
      for (int j=0; j<sos.rows(); j++)
        {

          double v0=0.0, v1=0.0, v2=0.0;

          double a0 = sos(j,3);
          double a1 = sos(j,4)/a0;
          double a2 = sos(j,5)/a0;
          double b0 = sos(j,0)/a0;
          double b1 = sos(j,1)/a0;
          double b2 = sos(j,2)/a0;

          for (int i=0; i<n; i++)
            {
              v0=x(i,k)-a1*v1-a2*v2;
              y(i,k)=b0*v0+b1*v1+b2*v2;
              v2=v1;
              v1=v0;
            }

          x.insert(y.column(k),0,k);

        }
    }

  if (isrowvector)
    y=y.transpose();

  retval(0)=y;

  return retval;
}
