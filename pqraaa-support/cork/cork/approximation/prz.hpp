//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_prz_hpp
#define cork_approximation_prz_hpp

#include <cork/approximation/prz_options.hpp>
#include <cork/approximation/prz_approximation.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>

namespace CORK { namespace approximation {

  template <typename T, typename R, typename Z, typename F, typename W>
  prz_approximation<T> PRZ( R const& rat, Z const& zj, F const& fj, W const& wj ) {
    function [pol, res, zer] = prz(r, zj, fj, wj)
      % Compute poles, residues, and zeros of rational function in barycentric form.
      m = length(wj);

    % Compute poles via generalized eigenvalue problem:
      B = eye(m+1);
    B(1,1) = 0;
    E = [0 wj.'; ones(m, 1) diag(zj)];
    pol = eig(E, B);
    % Remove zeros of denominator at infinity:
      pol = pol(~isinf(pol));

    % Compute residues via discretized Cauchy integral:
      dz = 1e-5*exp(2i*pi*(1:4)/4);

    pp = pol+dz;
    rvals = r(pp(:));
    res = zeros(length(pol),size(rvals,2));
    for it = 1:size(rvals,2)
          res(:,it) = reshape(rvals(:,it),[],4)*dz.'/4;
    end

      % Compute zeros via generalized eigenvalue problem:
      for it = 1:size(fj,2)
            E = [0 (wj.*fj(:,it)).'; ones(m, 1) diag(zj)];
        zer{it} = eig(E, B);
            % Remove zeros of numerator at infinity:
                  zer{it} = zer{it}(~isinf(zer{it}));
            end
              end % End of PRZ().
  } // PRZ()

} } // namespace CORK::approximation

#endif
