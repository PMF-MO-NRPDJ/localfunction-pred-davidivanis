#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

// Lagrangeovi konačni elementi
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include "interp_error.hh"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    // Grid
    const int dim = 2;
    using GridType = Dune::YaspGrid<dim>;
    using GridView = typename GridType::LeafGridView;
    Dune::FieldVector<double, dim> L(5);  // Domena = (0,5)x(0,2)
    L[1] = 2;
    std::array<int,dim> N{25,10};
    GridType grid(L,N);
    auto gv = grid.leafGridView();
    double h = 0.2;  // inicijalni korak mreže

    // Konačni elementi
    const int polDeg = 2;
    using FEM = Dune::LagrangeCubeLocalFiniteElement<double, double, dim, polDeg>;
    // za simplekse:
    //using FEM = Dune::LagrangeSimplexLocalFiniteElement<double, double, dim, polDeg>;
    FEM fem;

    // Struktura LocalFiniteElemet objekta. Podtipovi
    using LocalBasisType         = FEM::Traits::LocalBasisType;
    using LocalCoefficientsType  = FEM::Traits::LocalCoefficientsType;
    using LocalInterpolationType = FEM::Traits::LocalInterpolationType;

    // Struktura LocalFiniteElemet objekta. Podobjekti.
    LocalBasisType         const & basis = fem.localBasis();
    LocalCoefficientsType  const & coeff = fem.localCoefficients();
    LocalInterpolationType const & inter = fem.localInterpolation();

    // Lokalna baza
    // Tipovi koji su definirani u lokalnoj bazi
    using DomainType   = LocalBasisType::Traits::DomainType; // Tip varijable u domeni funkcije
    using RangeType    = LocalBasisType::Traits::RangeType;  // Tip varijable u kodomeni
    using JacobianType = LocalBasisType::Traits::JacobianType; // Tip jakobijana

    DomainType x; // DomainType je FieldVector<double,dim>
    std::vector<RangeType> phi;// RangeType je "ekvivalentan" s double jer su funkcije skalarne.
    std::vector<JacobianType> grad_phi;

    x[0]=0.5; x[1]=0.5;
    basis.evaluateFunction(x,phi);      // izračunaj sve bazne funkcije u točki x
    basis.evaluateJacobian(x,grad_phi); // izračunaj gradijente svih baznih funkcija u točki x

    //std::array<unsigned int, dim> order{2,1};
    //std::vector<RangeType> parc;
    //basis.partial(order,x,parc);

    std::cout << "=== Vrijenosti baznih funkcija i njenih derivacija ===\n";
    for(unsigned int i=0; i<phi.size(); ++i){
        std::cout << "phi_"<< i << "(" << x <<") = " << phi[i] << ", ";
        std::cout << "grad phi_"<< i << "(" << x << ") = " << "(" << grad_phi[i][0]<<")" << std::endl;

    //    std::cout << "D_{2,1}phi_"<< i << "(" << x <<") = " << parc[i] << "\n";
    }

    // Lokalni koeficijenti
    //coeff.localKey(i) -> [index, codim, k]
    std::cout << "== Lokalni koeficijenti ===\n";
    for(unsigned int n_deg=0; n_deg<coeff.size(); ++n_deg){
        auto key = coeff.localKey(n_deg);
        if     (key.codim() == dim)   std::cout << "vrh_";
        else if(key.codim() == dim-1) std::cout << "brid_";
        else if(key.codim() == 0)     std::cout << "element_";
        std::cout << key.subEntity(); // broj subentiteta
        std::cout << ", k= " << key.index();
        std::cout << "\n";
    }

    // Lokalna interpolacija
    std::cout << "== Lokalna interpolacija ===\n";
    std::vector<RangeType> coefficients;
    auto funkc = [](DomainType x){ return x*x;};
    inter.interpolate(funkc, coefficients);
    // Ispišimo keficijente
    std::cout << "Koeficijenti : ";
    for(double val : coefficients) std::cout << val <<", ";
    std::cout << std::endl;

    // Zadatak s greškom interpolacije.
    // Prevoditelj ne može deducirati Element.
    // using Element = GridView::template Codim<0>::Entity;
    double error_prev = interpolationError(gv, fem);
    for(int i=0; i<4; ++i){
        grid.globalRefine(1);
        double error = interpolationError(gv, fem);
        double alpha = std::log(error_prev/error) / std::log(2.0);
        double C = error_prev/std::pow(h,alpha);
        std::cout << "Error = " << error
                  << ", alpha = " << alpha << ", C = " << C << "\n";
        h /= 2;
        error_prev = error;
    }
    return 0;
}
