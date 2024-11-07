//  (C) Copyright Karl Meerbergen 2011. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)
//
//  This files uses gnuplot-gcc which is protected under GNU GPL-v2.
//  By using GLAS, the user accepts the terms of licence of gnuplot-gcc.

#ifndef glas2_gnuplot_gnuplot_hpp
#define glas2_gnuplot_gnuplot_hpp

#include <glas2/gnuplot/wait.hpp>
#include <glas2/gnuplot/impl/gnuplot_i.hpp>
#include <glas2/gnuplot/glas2std.hpp>
#include <glas2/matrix/container/shared_matrix.hpp>
#include <cassert>

namespace glas2 { namespace gnuplot {

class gnuplot
: public Gnuplot
{
    public:

	    ///\brief set a style during construction
        gnuplot(const std::string &style = "points")
        : Gnuplot( style )
        {}

        /// plot a single std::vector at one go
        template <typename X>
        gnuplot(const X &x,
                const std::string &title = "",
                const std::string &style = "points",
                const std::string &labelx = "x",
                const std::string &labely = "y")
        : Gnuplot( style )
        {
            set_xlabel(labelx);
            set_ylabel(labely);

            plot_x(x,title);
        }

         /// plot pairs std::vector at one go
        template <typename X, typename Y>
        gnuplot(const X &x, const Y &y,
                const std::string &title = "",
                const std::string &style = "points",
                const std::string &labelx = "x",
                const std::string &labely = "y")
        : Gnuplot( style )
        {
            set_xlabel(labelx);
            set_ylabel(labely);

            plot_xy(x,y,title);
        }

         /// plot triples std::vector at one go
        template <typename X, typename Y, typename Z>
        gnuplot(const X &x, const Y &y, const Z &z,
                const std::string &title = "",
                const std::string &style = "points",
                const std::string &labelx = "x",
                const std::string &labely = "y",
                const std::string &labelz = "z")
        : Gnuplot( style )
        {
            set_xlabel(labelx);
            set_ylabel(labely);
            set_zlabel(labelz);

            plot_xyz(x,y,z,title);
        }

    //----------------------------------------------------------------------------------

    /// send a command to gnuplot
    gnuplot& cmd(const std::string &cmdstr) { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).cmd( cmdstr ) ) ; }
	// ---------------------------------------------------------------------------------
	///\brief Sends a command to an active gnuplot session, identical to cmd()
	/// send a command to gnuplot using the <<  operator
	///
	/// \param cmdstr --> the command string
	/// 
	/// \return <-- a reference to the gnuplot object	
	// ---------------------------------------------------------------------------------
    inline gnuplot& operator<<(const std::string &cmdstr){
        cmd(cmdstr);
        return(*this);
    }



    //----------------------------------------------------------------------------------
    // show on screen or write to file

    /// sets terminal type to terminal_std
    gnuplot& showonscreen() { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).showonscreen() ) ; }

    /// saves a gnuplot session to a postscript file, filename without extension
    //Gnuplot& savetops(const std::string &filename = "gnuplot_output");
    gnuplot& print(const std::string &filename, const std::string fmt="png")
    {
      {
        std::ostringstream cmdstr;
        cmdstr << std::string("set terminal ") ;
        cmdstr << fmt;
        cmd(cmdstr.str());
      }

      {
        std::ostringstream cmdstr;
        cmdstr << std::string("set output \"") << filename << std::string("\"");
        cmd(cmdstr.str());
        cmd("replot");
      }

      return *this;
    }


    //----------------------------------------------------------------------------------
    // set and unset

    /// set line style (some of these styles require additional information):
    ///  lines, points, linespoints, impulses, dots, steps, fsteps, histeps,
    ///  boxes, histograms, filledcurves
    gnuplot& style(const std::string &stylestr = "points") { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).set_style( stylestr ) ) ; }

    /// interpolation and approximation of data, arguments:
    ///  csplines, bezier, acsplines (for data values > 0), sbezier, unique, frequency
    /// (works only with plot_x, plot_xy, plotfile_x, plotfile_xy
    /// (if smooth is set, set_style has no effekt on data plotting)
    //Gnuplot& set_smooth(const std::string &stylestr = "csplines");

    // ----------------------------------------------------------------------
    /// \brief unset smooth
    /// attention: smooth is not set by default
    /// 
    /// \param ---
    /// 
    /// \return <-- a reference to a gnuplot object
    // ----------------------------------------------------------------------
    //inline Gnuplot& unset_smooth(){ smooth = ""; return *this;}; 


    /// scales the size of the points used in plots
    gnuplot& set_pointsize(const double pointsize = 1.0) { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).set_pointsize(pointsize) ) ; }

    /// turns grid on/off
    //inline Gnuplot& set_grid()	{cmd("set grid");return *this;};
    /// grid is not set by default	
    //inline Gnuplot& unset_grid(){cmd("unset grid");return *this;}; 

    // -----------------------------------------------
    /// set the mulitplot mode
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------
    //inline Gnuplot& set_multiplot(){cmd("set multiplot") ;return *this;};

    // -----------------------------------------------
    /// unsets the mulitplot mode
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------
    //inline Gnuplot& unset_multiplot(){cmd("unset multiplot");return *this;};
    


    /// set sampling rate of functions, or for interpolating data
    //Gnuplot& set_samples(const int samples = 100);
    /// set isoline density (grid) for plotting functions as surfaces (for 3d plots)
    //Gnuplot& set_isosamples(const int isolines = 10);

    // --------------------------------------------------------------------------
    /// enables/disables hidden line removal for surface plotting (for 3d plot)
    ///
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // --------------------------------------------------------------------------
    //Gnuplot& set_hidden3d(){cmd("set hidden3d");return *this;};

    // ---------------------------------------------------------------------------
    /// hidden3d is not set by default
    ///
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // ---------------------------------------------------------------------------
    //inline Gnuplot& unset_hidden3d(){cmd("unset hidden3d"); return *this;}; 

    /// enables/disables contour drawing for surfaces (for 3d plot)
    ///  base, surface, both
    //Gnuplot& set_contour(const std::string &position = "base");
    // --------------------------------------------------------------------------
    /// contour is not set by default, it disables contour drawing for surfaces
    ///
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // ------------------------------------------------------------------
    //inline Gnuplot& unset_contour(){cmd("unset contour");return *this;};

    // ------------------------------------------------------------
    /// enables/disables the display of surfaces (for 3d plot)
    ///
    /// \param ---	
    ///
    /// \return <-- reference to the gnuplot object
    // ------------------------------------------------------------------
    //inline Gnuplot& set_surface(){cmd("set surface");return *this;}; 

    // ----------------------------------------------------------
    /// surface is set by default,
    /// it disables the display of surfaces (for 3d plot)
    ///
    /// \param ---	
    ///
    /// \return <-- reference to the gnuplot object
    // ------------------------------------------------------------------
    //inline Gnuplot& unset_surface(){cmd("unset surface"); return *this;}


    /// switches legend on/off
    /// position: inside/outside, left/center/right, top/center/bottom, nobox/box
    //Gnuplot& set_legend(const std::string &position = "default"); 

    // ------------------------------------------------------------------
    /// \brief  Switches legend off
    /// attention:legend is set by default
    ///
    /// \param ---	
    ///
    /// \return <-- reference to the gnuplot object
    // ------------------------------------------------------------------
    //inline Gnuplot& unset_legend(){cmd("unset key"); return *this;}

    // ----------------------------------------------------------------------- 
    /// \brief sets and clears the title of a gnuplot session
    ///
    /// \param title --> the title of the plot [optional, default == ""]
    ///
    /// \return <-- reference to the gnuplot object
    // ----------------------------------------------------------------------- 
    inline gnuplot& title(const std::string &title = "") { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).set_title( title ) ) ; }

    //----------------------------------------------------------------------------------
    ///\brief Clears the title of a gnuplot session
    /// The title is not set by default.
    ///
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // ---------------------------------------------------------------------------------
    inline gnuplot& unset_title() {
      static_cast<Gnuplot&>(*this).unset_title() ;
      return *this ;
    }


    /// set x axis label
    //Gnuplot& set_ylabel(const std::string &label = "x");
    /// set y axis label
    //Gnuplot& set_xlabel(const std::string &label = "y");
    /// set z axis label
    //Gnuplot& set_zlabel(const std::string &label = "z");

    /// set axis - ranges
    /*Gnuplot& set_xrange(const double iFrom,
                        const double iTo);*/
    /// set y-axis - ranges
    /*Gnuplot& set_yrange(const double iFrom,
                        const double iTo);*/
    /// set z-axis - ranges
    /*Gnuplot& set_zrange(const double iFrom,
                        const double iTo);*/
    /// autoscale axis (set by default) of xaxis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------
    //inline Gnuplot& set_xautoscale(){cmd("set xrange restore");cmd("set autoscale x");return *this;};

    // -----------------------------------------------
    /// autoscale axis (set by default) of yaxis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------
    //inline Gnuplot& set_yautoscale(){cmd("set yrange restore");cmd("set autoscale y");return *this;};

    // -----------------------------------------------
    /// autoscale axis (set by default) of zaxis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------
    //inline Gnuplot& set_zautoscale(){cmd("set zrange restore");cmd("set autoscale z");return *this;};


    /// turns on/off log scaling for the specified xaxis (logscale is not set by default)
    //Gnuplot& set_xlogscale(const double base = 10);
    /// turns on/off log scaling for the specified yaxis (logscale is not set by default)
    //Gnuplot& set_ylogscale(const double base = 10);
    /// turns on/off log scaling for the specified zaxis (logscale is not set by default)
    //Gnuplot& set_zlogscale(const double base = 10);

    // ----------------------------------------------- 
    /// turns off log scaling for the x axis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------	
    //inline Gnuplot& unset_xlogscale(){cmd("unset logscale x"); return *this;};

    // ----------------------------------------------- 
    /// turns off log scaling for the y axis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------	
    //inline Gnuplot& unset_ylogscale(){cmd("unset logscale y"); return *this;};

    // ----------------------------------------------- 
    /// turns off log scaling for the z axis
    /// 
    /// \param ---
    ///
    /// \return <-- reference to the gnuplot object
    // -----------------------------------------------		
    //inline Gnuplot& unset_zlogscale(){cmd("unset logscale z"); return *this;};


    /// set palette range (autoscale by default)
    //Gnuplot& set_cbrange(const double iFrom, const double iTo);


    //----------------------------------------------------------------------------------
    // plot

    /// plot a single std::vector: x
    ///   from file
    /*Gnuplot& plotfile_x(const std::string &filename,
                        const unsigned int column = 1,
                        const std::string &title = "");*/
    ///   from std::vector
    template<typename X>
    gnuplot& plot_x(const X& x, const std::string &title = "")
    {
      static_cast<Gnuplot&>(*this).plot_x( glas2std(x), title ) ;
      return *this;
    }



    /// plot x,y pairs: x y
    ///   from file
    /*Gnuplot& plotfile_xy(const std::string &filename,
                         const unsigned int column_x = 1,
                         const unsigned int column_y = 2,
                         const std::string &title = "");*/
    ///   from data
    template<typename X, typename Y>
    gnuplot& plot_xy(const X& x, const Y& y, const std::string &title = ""){
      assert( x.size() == y.size() ) ;
      static_cast<Gnuplot&>(*this).plot_xy( glas2std(x), glas2std(y), title ) ;
      return *this ;
    }


    /// plot x,y pairs with dy errorbars: x y dy
    ///   from file
    /*Gnuplot& plotfile_xy_err(const std::string &filename,
                             const unsigned int column_x  = 1,
                             const unsigned int column_y  = 2,
                             const unsigned int column_dy = 3,
                             const std::string &title = "");
    ///   from data
    template<typename X, typename Y, typename E>
    Gnuplot& plot_xy_err(const X &x, const Y &y, const E &dy,
                         const std::string &title = "");*/


    /// plot x,y,z triples: x y z
    ///   from file
    /*gnuplot& plotfile_xyz(const std::string &filename,
                          const unsigned int column_x = 1,
                          const unsigned int column_y = 2,
                          const unsigned int column_z = 3,
                          const std::string &title = "");*/
    ///   from std::vector
    template<typename X, typename Y, typename Z>
    Gnuplot& plot_xyz(const X &x,
                      const Y &y,
                      const Z &z,
                      const std::string &title = "")
    {
      assert( x.size() == y.size() ) ;
      assert( x.size() == z.size() ) ;
      static_cast<Gnuplot&>(*this).plot_xyz( glas2std(x), glas2std(y), glas2std(z), title ) ;
      return *this ;
    }



    /// plot an equation of the form: y = ax + b, you supply a and b
   /* Gnuplot& plot_slope(const double a,
                        const double b,
                        const std::string &title = "");*/


    /// plot an equation supplied as a std::string y=f(x), write only the function f(x) not y=
    /// the independent variable has to be x
    /// binary operators: ** exponentiation, * multiply, / divide, + add, - substract, % modulo
    /// unary operators: - minus, ! factorial
    /// elementary functions: rand(x), abs(x), sgn(x), ceil(x), floor(x), int(x), imag(x), real(x), arg(x),
    ///   sqrt(x), exp(x), log(x), log10(x), sin(x), cos(x), tan(x), asin(x), acos(x), atan(x), atan2(y,x),
    ///   sinh(x), cosh(x), tanh(x), asinh(x), acosh(x), atanh(x)
    /// special functions: erf(x), erfc(x), inverf(x), gamma(x), igamma(a,x), lgamma(x), ibeta(p,q,x),
    ///   besj0(x), besj1(x), besy0(x), besy1(x), lambertw(x)
    /// statistical fuctions: norm(x), invnorm(x)
    /*
    Gnuplot& plot_equation(const std::string &equation,
                           const std::string &title = "");

    /// plot an equation supplied as a std::string z=f(x,y), write only the function f(x,y) not z=
    /// the independent variables have to be x and y
    Gnuplot& plot_equation3d(const std::string &equation,
                             const std::string &title = "");


    /// plot image
    Gnuplot& plot_image(const unsigned char *ucPicBuf,
                        const unsigned int iWidth,
                        const unsigned int iHeight,
                        const std::string &title = "");
*/

    //----------------------------------------------------------------------------------
    ///\brief replot repeats the last plot or splot command.
    ///  this can be useful for viewing a plot with different set options,
    ///  or when generating the same plot for several devices (showonscreen, savetops)
    /// 
    /// \param ---
    ///
    /// \return ---
    //----------------------------------------------------------------------------------
    inline gnuplot& replot(void){
      static_cast<Gnuplot&>(*this).replot() ;
      return *this ;
    }

    /// resets a gnuplot session (next plot will erase previous ones)
    gnuplot& reset_plot() { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).reset_plot() ) ; }

    /// resets a gnuplot session and sets all variables to default
    gnuplot& reset_all() { return reinterpret_cast<gnuplot&>( static_cast<Gnuplot&>(*this).reset_all() ) ; }

    template <typename X, typename Y, typename Z>
    gnuplot& surface( X const& x, Y const& y, Z const& z ) {
      assert( x.size()==z.num_rows() ) ;
      assert( y.size()==z.num_columns() ) ;
      std::ofstream tmp;
      std::string filename = create_tmpfile(tmp);
      if (filename == "")
        return *this;

      std::ostringstream cmdstr;
      cmdstr << std::string("splot \"") << filename << std::string("\" with ") << pstyle;

      glas2::shared_matrix<typename Z::value_type> z_temp( z.num_rows(), z.num_columns() ) ;
      z_temp = z ;

      for ( int i = 0; i < x.size(); i++) {
        for ( int j = 0; j < y.size(); j++) {
          tmp << x(i) << std::string(" ") << y(j) << std::string(" ") << z_temp(i,j) <<std::endl;
        }
        tmp << std::endl;
      }

      this->cmd(cmdstr.str());

      return *this ;
    }

};

} } // namespace glas2::gnuplot

#endif
