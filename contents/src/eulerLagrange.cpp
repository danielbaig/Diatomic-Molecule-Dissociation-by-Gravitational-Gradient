// Configured==true
#include "eulerLagrange.h"

const double epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


std::tuple<double, double, double, double> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const double sigma_p, const double delta_p, const double chi_p, const double dlambda)
{
	double t_new{};
	double r_new{};
	double phi_new{};
	double theta_new{};

	const double tdot{ (p->t - p->t_prev) / dlambda };
	const double rdot{ (p->r - p->r_prev) / dlambda };
	const double phidot{ (p->phi - p->phi_prev) / dlambda };
	const double thetadot{ (p->theta - p->theta_prev) / dlambda };

	const double tsdot{ (ps->t - ps->t_prev) / dlambda };
	const double rsdot{ (ps->r - ps->r_prev) / dlambda };
	const double phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const double thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const double tsddot{(ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const double rsddot{(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const double phisddot{(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const double thetasddot{(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};
	const double nx{p->r*sin(p->theta)*cos(p->phi) - ps->r*sin(ps->theta)*cos(ps->phi)};
	const double ny{p->r*sin(p->theta)*sin(p->phi) - ps->r*sin(ps->theta)*sin(ps->phi)};
	const double nz{p->r*cos(p->theta) - ps->r*cos(ps->theta)};
	const double nmag{sqrt(nx*nx + ny*ny + nz*nz)};

	t_new = dlambda*(phidot*chi_p + (p->energy*sigma_p)/(2.*delta_p));
	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((0.5*delta_p*(2.*p->r*(pow(thetadot,2.) + pow(rdot,2.)/delta_p) + (2.*p->r*pow(tdot - phidot*chi_p,2.)*delta_p)/pow(sigma_p,2.) + (2.*(BH->mass - p->r)*pow(tdot - phidot*chi_p,2.))/sigma_p - (2.*(BH->mass - p->r)*pow(rdot,2.)*sigma_p)/pow(delta_p,2.) - (2.*pow(phidot,2.)*pow(p->r,5)*pow(sin(p->theta),2.))/pow(sigma_p,2.) - (2.*rdot*(2.*p->r*rdot - 2.*thetadot*(BH->l + BH->a*cos(p->theta))*sin(p->theta)))/delta_p))/sigma_p);
	phi_new = dlambda*((2.*tdot*chi_p*delta_p - p->L*sigma_p)/(2.*pow(chi_p,2.)*delta_p - 2.*pow(p->r,4)*pow(sin(p->theta),2.)));
	theta_new = 2 * p->theta - p->theta_prev + dlambda*dlambda*((0.5*(-4.*p->r*rdot*thetadot + 4.*pow(thetadot,2.)*(BH->l + BH->a*cos(p->theta))*sin(p->theta) - 2.*(BH->l + BH->a*cos(p->theta))*(pow(thetadot,2.) + pow(rdot,2.)/delta_p)*sin(p->theta) - (2.*pow(tdot - phidot*chi_p,2.)*(BH->l + BH->a*cos(p->theta))*delta_p*sin(p->theta))/pow(sigma_p,2.) + (2.*pow(phidot,2.)*pow(p->r,4)*(BH->l + BH->a*cos(p->theta))*pow(sin(p->theta),3))/pow(sigma_p,2.) + (pow(phidot,2.)*pow(p->r,4)*sin(2.*p->theta))/sigma_p + (2.*phidot*(tdot - phidot*chi_p)*delta_p*(2.*BH->l*sin(p->theta) + BH->a*sin(2.*p->theta)))/sigma_p))/sigma_p);

	return {t_new, r_new, phi_new, theta_new};
}
