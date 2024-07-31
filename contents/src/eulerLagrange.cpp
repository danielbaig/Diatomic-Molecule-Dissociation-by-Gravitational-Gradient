// Configured==true
#include "eulerLagrange.h"

const cpp_dec_float_25 epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


std::tuple<cpp_dec_float_25, cpp_dec_float_25,	cpp_dec_float_25, cpp_dec_float_25> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 t_new{};
	cpp_dec_float_25 r_new{};
	cpp_dec_float_25 phi_new{};
	cpp_dec_float_25 theta_new{};

	const cpp_dec_float_25 sigma_p2{sigma_p * sigma_p};
	const cpp_dec_float_25 sigma_p3{sigma_p * sigma_p2};
	const cpp_dec_float_25 sigma_p4{sigma_p2 * sigma_p2};
	const cpp_dec_float_25 delta_p2{delta_p * delta_p};
	const cpp_dec_float_25 delta_p3{delta_p * delta_p2};
	const cpp_dec_float_25 delta_p4{delta_p2 * delta_p2};
	const cpp_dec_float_25 chi_p2{chi_p * chi_p};
	const cpp_dec_float_25 chi_p3{chi_p * chi_p2};
	const cpp_dec_float_25 chi_p4{chi_p2 * chi_p2};

	const cpp_dec_float_25 tdot{ (p->t - p->t_prev) / dlambda };
	const cpp_dec_float_25 rdot{ (p->r - p->r_prev) / dlambda };
	const cpp_dec_float_25 phidot{ (p->phi - p->phi_prev) / dlambda };
	const cpp_dec_float_25 thetadot{ (p->theta - p->theta_prev) / dlambda };

	const cpp_dec_float_25 tsdot{ (ps->t - ps->t_prev) / dlambda };
	const cpp_dec_float_25 rsdot{ (ps->r - ps->r_prev) / dlambda };
	const cpp_dec_float_25 phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const cpp_dec_float_25 thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const cpp_dec_float_25 tsddot{(ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 rsddot{(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 phisddot{(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 thetasddot{	(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};

	const cpp_dec_float_25 nx{p->r*sin(p->theta)*cos(p->phi) - ps->r*sin(ps->theta)*cos(ps->phi)};
	const cpp_dec_float_25 ny{p->r*sin(p->theta)*sin(p->phi) - ps->r*sin(ps->theta)*sin(ps->phi)};
	const cpp_dec_float_25 nz{p->r*cos(p->theta) - ps->r*cos(ps->theta)};
	const cpp_dec_float_25 nmag{sqrt(nx*nx + ny*ny + nz*nz)};
	const cpp_dec_float_25 nmag2{nmag*nmag};
	const cpp_dec_float_25 nmag5{nmag*nmag2*nmag2};
	const cpp_dec_float_25 nmag6{nmag5*nmag};
	const cpp_dec_float_25 nmag7{nmag6*nmag};
	const cpp_dec_float_25 nmag8{nmag6*nmag2};
	const cpp_dec_float_25 nmag12{nmag6*nmag6};
	const cpp_dec_float_25 nmag13{nmag6*nmag7};
	const cpp_dec_float_25 nmag14{nmag7*nmag7};

	const cpp_dec_float_25 ral_sqr{BH->a2 + BH->l2 + p->r2};
	const cpp_dec_float_25 sintheta{sin(p->theta)};
	const cpp_dec_float_25 sintheta2{sintheta*sintheta};
	const cpp_dec_float_25 sintheta3{sintheta2*sintheta};
	const cpp_dec_float_25 sintheta4{sintheta2*sintheta2};
	const cpp_dec_float_25 sintheta5{sintheta2*sintheta3};
	const cpp_dec_float_25 sintheta6{sintheta3*sintheta3};
	const cpp_dec_float_25 sintheta7{sintheta3*sintheta4};

	const cpp_dec_float_25 costheta{cos(p->theta)};
	const cpp_dec_float_25 sin2theta{sin(2.*p->theta)};
	const cpp_dec_float_25 sinphi{sin(p->phi)};
	const cpp_dec_float_25 cosphi{cos(p->phi)};

	const cpp_dec_float_25 sinthetas{sin(ps->theta)};
	const cpp_dec_float_25 costhetas{cos(ps->theta)};
	const cpp_dec_float_25 sinphis{sin(ps->phi)};
	const cpp_dec_float_25 cosphis{cos(ps->phi)};
	const cpp_dec_float_25 sinthetasdot{sin(thetasdot)};
	const cpp_dec_float_25 costhetasdot{cos(thetasdot)};
	const cpp_dec_float_25 sinphisdot{sin(phisdot)};
	const cpp_dec_float_25 cosphisdot{cos(phisdot)};
	const cpp_dec_float_25 sinthetas2{sinthetas*sinthetas};
	const cpp_dec_float_25 costhetas2{costhetas*costhetas};
	const cpp_dec_float_25 sinphis2{sinphis*sinphis};
	const cpp_dec_float_25 cosphis2{cosphis*cosphis};
	const cpp_dec_float_25 rsdot2{rsdot*rsdot};
	const cpp_dec_float_25 phisdot2{phisdot*phisdot};

	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-0.5*((-24.*ps->epsilon*ps->sigma6*((thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag + (rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag + (phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag + (thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag + (rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag + (phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag7) - (48.*ps->epsilon*ps->sigma6*(-nmag6 + ps->sigma6)*((thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag + (rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag + (phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag + (thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag + (rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag + (phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13)));
	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*(0.5*(2.*p->r*((thetadot*thetadot) + (phidot*phidot)*sintheta2) + (8.*ps->epsilon*p->r*ps->sigma6*(-nmag6 + ps->sigma6)*(thetadot*thetasdot + phidot*phisdot*sintheta2))/(p->mass*nmag12) - (24.*ps->epsilon*ps->sigma6*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta)*(rdot*rsdot - tdot + p->r2*(thetadot*thetasdot + phidot*phisdot*sintheta2)))/(p->mass*nmag8) - (48.*ps->epsilon*ps->sigma6*(-nmag6 + ps->sigma6)*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta)*(rdot*rsdot - tdot + p->r2*(thetadot*thetasdot + phidot*phisdot*sintheta2)))/(p->mass*nmag14) - (ps->epsilon*ps->sigma6*(-4.*nmag7*rsddot + 4.*nmag*rsddot*ps->sigma6 + rsdot*ps->sigma6*((-48.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag - (48.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag - (48.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag - (48.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag - (48.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag - (48.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag) + nmag6*rsdot*((24.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag + (24.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag + (24.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag + (24.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag + (24.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag + (24.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag)))/(p->mass*nmag13)));
	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((0.5*pow(1. / sin(p->theta),2)*(-(phidot*p->r*sintheta*(4.*p->r*thetadot*costheta + 4.*rdot*sintheta)) - (ps->epsilon*p->r*ps->sigma6*sintheta*(-8.*phisdot*rdot*sintheta + p->r*(-8.*phisdot*thetadot*costheta - 4.*phisddot*sintheta)))/(p->mass*nmag6) - (ps->epsilon*p->r*ps->sigma12*sintheta*(8.*phisdot*rdot*sintheta + p->r*(8.*phisdot*thetadot*costheta + 4.*phisddot*sintheta)))/(p->mass*nmag12) - (48.*ps->epsilon*ps->sigma6*(-0.5*nmag6 + ps->sigma6)*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta)*(rdot*rsdot - tdot + p->r2*(thetadot*thetasdot + phidot*phisdot*sintheta2)))/(p->mass*nmag14) - (ps->epsilon*phisdot*p->r2*ps->sigma12*sintheta2*((-48.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag - (48.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag - (48.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag - (48.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag - (48.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag - (48.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13) - (ps->epsilon*phisdot*p->r2*ps->sigma6*sintheta2*((24.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag + (24.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag + (24.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag + (24.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag + (24.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag + (24.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag7)))/p->r2);
	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*((0.5*(-4.*p->r*rdot*thetadot - (ps->epsilon*p->r*ps->sigma6*(-4.*p->r*thetasddot - 8.*rdot*thetasdot))/(p->mass*nmag6) - (ps->epsilon*p->r*ps->sigma12*(4.*p->r*thetasddot + 8.*rdot*thetasdot))/(p->mass*nmag12) - (24.*ps->epsilon*ps->sigma6*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta)*(rdot*rsdot - tdot + p->r2*(thetadot*thetasdot + phidot*phisdot*sintheta2)))/(p->mass*nmag8) - (48.*ps->epsilon*ps->sigma6*(-nmag6 + ps->sigma6)*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta)*(rdot*rsdot - tdot + p->r2*(thetadot*thetasdot + phidot*phisdot*sintheta2)))/(p->mass*nmag14) + (phidot*phidot)*p->r2*sin2theta + (4.*ps->epsilon*phidot*phisdot*p->r2*ps->sigma6*(-nmag6 + ps->sigma6)*sin2theta)/(p->mass*nmag12) - (ps->epsilon*p->r2*ps->sigma12*thetasdot*((-48.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag - (48.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag - (48.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag - (48.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag - (48.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag - (48.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13) - (ps->epsilon*p->r2*ps->sigma6*thetasdot*((24.*thetadot*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/nmag + (24.*rdot*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta))/nmag + (24.*phidot*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta))/nmag + (24.*thetasdot*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas))/nmag + (24.*rsdot*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas))/nmag + (24.*phisdot*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag7)))/p->r2);

	return {t_new, r_new, phi_new, theta_new};
}
