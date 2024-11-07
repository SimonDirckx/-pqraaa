#ifndef POINTMAT_H
#define POINTMAT_H

#include <grid/grid.hh>
#include <Eigen/Dense>
#include <complex>

namespace Grid
{
    class pointMat
    {
    private:
        const Grid& _grid;

    public:
        pointMat(const Grid& grid) : _grid(grid) {}

        Eigen::MatrixXcd computeKernelMat(double kappa)
        {
            // Get the vertices of the grid
            const auto& vertices = _grid.getVertices();
            
            // Initialize the kernel matrix
            Eigen::MatrixXcd kernelMat(vertices.cols(), vertices.cols());

            // Compute the kernel matrix
            for (int i = 0; i < vertices.cols(); ++i)
            {
                for (int j = 0; j < vertices.cols(); ++j)
                {
                    double distance = (vertices.col(i) - vertices.col(j)).norm();
                    std::complex<double> value = std::exp(std::complex<double>(0, kappa * distance));
                    kernelMat(i, j) = value;
                }
            }

            return kernelMat;
        }
    };
}

#endif // POINTMAT_H
