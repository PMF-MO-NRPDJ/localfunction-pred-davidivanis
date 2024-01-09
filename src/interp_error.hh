#pragma once
#include <dune/geometry/quadraturerules.hh>

template <int dim>
using Point = Dune::FieldVector<double, dim>;

// Funkcija čiji interpoland računamo zadana u globalnim koordinatama.
template <int dim>
double funkcija(Point<dim> const & x_global)
{
    double rez = 1.0;
    for(int i=0; i<dim; ++i)
        rez *= std::sin(x_global[i]);
    return rez;
}

// Lokalna verzija gornje funkcije. Uzima element u konstruktoru
// i lokalnu koordinatu u operatoru funkcijskog poziva
// (umjesto globalne koordinate).
template <typename Element>
struct funkcija_loc
{
    static const int dim = Element::dimension;
    Element element;

    funkcija_loc(Element elem) : element{elem} {}

    double operator()(Point<dim> const & x_local) const{
        auto x_global = element.geometry().global(x_local);
        return funkcija(x_global);
    }
};



template <typename GridView, typename FEM>
double interpolationError(GridView const & gv, FEM const & fem)
{
    using RangeType = typename FEM::Traits::LocalBasisType::Traits::RangeType;
    const int dim = GridView::dimension;
    auto const & basis = fem.localBasis();
    auto const & interp = fem.localInterpolation();
    std::vector<RangeType> phi;
    std::vector<double> coeff;

    double global_l2_error2 = 0.0;
    for(auto const & elem : elements(gv))
    {
        // NASTAVI ....
        // formula numeričke integracije
        auto const & rule = Dune::QuadratureRules<double,dim>::rule(elem.geometry().type(), 6);
        double local_l2_error2 = 0.0;

        for(auto const & qpoint : rule)
        {
            const double wi = qpoint.weight();
            auto const & xi = qpoint.position();  // lokalna koordinata
            const double dx = elem.geometry().integrationElement(xi);
           // NASTAVI ...
        }
        global_l2_error2 += local_l2_error2;
    }
    return std::sqrt(global_l2_error2);
}
