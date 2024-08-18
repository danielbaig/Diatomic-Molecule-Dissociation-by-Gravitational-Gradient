// Configured==true
#include "eulerLagrange.h"

const cpp_dec_float_25 epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


cpp_dec_float_25 calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 t_new{};
	const cpp_dec_float_25 sigma_p2{ sigma_p*sigma_p };
	const cpp_dec_float_25 sigma_p3{ sigma_p*sigma_p2 };
	const cpp_dec_float_25 sigma_p4{ sigma_p2*sigma_p2 };
	const cpp_dec_float_25 delta_p2{ delta_p*delta_p };
	const cpp_dec_float_25 delta_p3{ delta_p*delta_p2 };
	const cpp_dec_float_25 delta_p4{ delta_p2*delta_p2 };
	const cpp_dec_float_25 chi_p2{ chi_p*chi_p };
	const cpp_dec_float_25 chi_p3{ chi_p*chi_p2 };
	const cpp_dec_float_25 chi_p4{ chi_p2*chi_p2 };

	const cpp_dec_float_25 tdot{ (p->t - p->t_prev) / dlambda };
	const cpp_dec_float_25 rdot{ (p->r - p->r_prev) / dlambda };
	const cpp_dec_float_25 phidot{ (p->phi - p->phi_prev) / dlambda };
	const cpp_dec_float_25 thetadot{ (p->theta - p->theta_prev) / dlambda };

	const cpp_dec_float_25 tsdot{ (ps->t - ps->t_prev) / dlambda };
	const cpp_dec_float_25 rsdot{ (ps->r - ps->r_prev) / dlambda };
	const cpp_dec_float_25 phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const cpp_dec_float_25 thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const cpp_dec_float_25 tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};

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

	const cpp_dec_float_25 temp00{ ps->epsilon*p->r*rdot*ps->sigma6*(-nmag6 + ps->sigma6) };
	const cpp_dec_float_25 temp01{ ps->epsilon*(-2.*BH->mass + 2.*p->r)*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp02{ (-2.*BH->mass + 2.*p->r)*rdot };
	const cpp_dec_float_25 temp03{ ps->epsilon*ps->sigma6 };
	const cpp_dec_float_25 temp04{ costheta + nx*cosphi*sintheta + ny*sinphi };
	const cpp_dec_float_25 temp05{ sintheta - nx*p->r*sinphi };
	const cpp_dec_float_25 temp06{ (nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis  };
	const cpp_dec_float_25 temp07{  nx*cosphis*sinthetas  };
	const cpp_dec_float_25 temp08{ (ny*ps->r*cosphis*sinthetas)  };
	const cpp_dec_float_25 temp09{ temp03*(-nmag6 + ps->sigma6) };
	const cpp_dec_float_25 temp010{  (phidot*(ny*p->r*cosphi*temp05*sintheta))/nmag + (thetasdot*(-temp06+ nz*ps->r*sinthetas))/nmag + (rsdot*(-(nz*costhetas) -temp07- ny*sinphis*sinthetas))/nmag  };
	const cpp_dec_float_25 temp011{ (BH->a*(ral_sqr) - chi_p*delta_p) };


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-((-((2.*pow(ral_sqr,2))/sigma_p - (2.*chi_p2*delta_p)/sigma_p)*((-4.*BH->a*p->r*rdot*(-phidot*(ral_sqr) + BH->a*tdot))/sigma_p2 - (8.*temp00*(BH->a2 - delta_p))/(p->mass*nmag12*sigma_p2) + (4.*p->r*rdot*(tdot - phidot*chi_p)*delta_p)/sigma_p2 - (4.*temp01*(-nmag6 + ps->sigma6))/(p->mass*nmag12*sigma_p) + (2.*temp02*(-tdot + phidot*chi_p))/sigma_p - (24.*temp03*(BH->a2 - delta_p)*((rdot*(nz*temp04*sintheta))/nmag +temp010+ (phisdot*(-temp08+ nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag7*sigma_p) - (48.*temp09*(BH->a2 - delta_p)*((rdot*(nz*temp04*sintheta))/nmag +temp010+ (phisdot*(-temp08+ nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13*sigma_p)) + ((-2.*BH->a*(ral_sqr))/sigma_p + (2.*chi_p*delta_p)/sigma_p)*((-4.*p->r*(ral_sqr)*rdot*(phidot*(ral_sqr) - BH->a*tdot))/sigma_p2 + (4.*p->r*rdot*chi_p*(-tdot + phidot*chi_p)*delta_p)/sigma_p2 + (16.*temp00*(BH->a*(ral_sqr) - chi_p*delta_p))/(p->mass*nmag12*sigma_p2) + (8.*temp01*(-nmag6 + ps->sigma6)*chi_p)/(p->mass*nmag12*sigma_p) + (2.*temp02*chi_p*(tdot - phidot*chi_p))/sigma_p - (96.*temp03*(-0.5*nmag6 + ps->sigma6)*(tdot*(-0.5*BH->a2 + 0.5*delta_p) + phidot*(BH->a3 + BH->a*BH->l2 + BH->a*p->r2 - chi_p*delta_p))*(ny*p->r*cosphi*temp05*sintheta))/(p->mass*nmag14*sigma_p) + (48.*temp03*temp011*((rdot*(nz*temp04*sintheta))/nmag +temp010+ (phisdot*(-temp08+ nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag7*sigma_p) + (96.*temp09*temp011*((rdot*(nz*temp04*sintheta))/nmag +temp010+ (phisdot*(-temp08+ nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13*sigma_p)))/(pow((-2.*BH->a*(ral_sqr))/sigma_p + (2.*chi_p*delta_p)/sigma_p,2) - ((2.*BH->a2)/sigma_p - (2.*delta_p)/sigma_p)*((2.*pow(ral_sqr,2))/sigma_p - (2.*chi_p2*delta_p)/sigma_p))));
	return t_new;
}

cpp_dec_float_25 calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 r_new{};
	const cpp_dec_float_25 sigma_p2{ sigma_p*sigma_p };
	const cpp_dec_float_25 sigma_p3{ sigma_p*sigma_p2 };
	const cpp_dec_float_25 sigma_p4{ sigma_p2*sigma_p2 };
	const cpp_dec_float_25 delta_p2{ delta_p*delta_p };
	const cpp_dec_float_25 delta_p3{ delta_p*delta_p2 };
	const cpp_dec_float_25 delta_p4{ delta_p2*delta_p2 };
	const cpp_dec_float_25 chi_p2{ chi_p*chi_p };
	const cpp_dec_float_25 chi_p3{ chi_p*chi_p2 };
	const cpp_dec_float_25 chi_p4{ chi_p2*chi_p2 };

	const cpp_dec_float_25 tdot{ (p->t - p->t_prev) / dlambda };
	const cpp_dec_float_25 rdot{ (p->r - p->r_prev) / dlambda };
	const cpp_dec_float_25 phidot{ (p->phi - p->phi_prev) / dlambda };
	const cpp_dec_float_25 thetadot{ (p->theta - p->theta_prev) / dlambda };

	const cpp_dec_float_25 tsdot{ (ps->t - ps->t_prev) / dlambda };
	const cpp_dec_float_25 rsdot{ (ps->r - ps->r_prev) / dlambda };
	const cpp_dec_float_25 phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const cpp_dec_float_25 thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const cpp_dec_float_25 tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};

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

	const cpp_dec_float_25 temp10{ ps->sigma6*(-nmag6 + ps->sigma6) };
	const cpp_dec_float_25 temp11{ (tdot*(BH->a2 - delta_p) - 2.*phidot*(BH->a*(ral_sqr) - chi_p*delta_p)) };
	const cpp_dec_float_25 temp12{ costheta + nx*cosphi*sintheta + ny*sinphi };


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((-0.5*delta_p*((2.*p->r*(rdot*rdot))/delta_p + (2.*p->r*pow(phidot*(ral_sqr) - BH->a*tdot,2))/sigma_p2 - (2.*p->r*pow(tdot - phidot*chi_p,2)*delta_p)/sigma_p2 + (8.*ps->epsilon*p->r*temp10*(tdot*(BH->a2 - delta_p) - 2.*phidot*(BH->a*(ral_sqr) - chi_p*delta_p)))/(p->mass*nmag12*sigma_p2) + ((-2.*BH->mass + 2.*p->r)*pow(tdot - phidot*chi_p,2))/sigma_p - (8.*ps->epsilon*(-2.*BH->mass + 2.*p->r)*temp10*(-0.5*tdot + phidot*chi_p))/(p->mass*nmag12*sigma_p) - ((-2.*BH->mass + 2.*p->r)*(rdot*rdot)*sigma_p)/delta_p2 + (24.*ps->epsilon*ps->sigma6*temp11*(nz*temp12*sintheta))/(p->mass*nmag8*sigma_p) + (48.*ps->epsilon*temp10*temp11*(nz*temp12*sintheta))/(p->mass*nmag14*sigma_p)))/sigma_p);
	return r_new;
}

cpp_dec_float_25 calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 phi_new{};
	const cpp_dec_float_25 sigma_p2{ sigma_p*sigma_p };
	const cpp_dec_float_25 sigma_p3{ sigma_p*sigma_p2 };
	const cpp_dec_float_25 sigma_p4{ sigma_p2*sigma_p2 };
	const cpp_dec_float_25 delta_p2{ delta_p*delta_p };
	const cpp_dec_float_25 delta_p3{ delta_p*delta_p2 };
	const cpp_dec_float_25 delta_p4{ delta_p2*delta_p2 };
	const cpp_dec_float_25 chi_p2{ chi_p*chi_p };
	const cpp_dec_float_25 chi_p3{ chi_p*chi_p2 };
	const cpp_dec_float_25 chi_p4{ chi_p2*chi_p2 };

	const cpp_dec_float_25 tdot{ (p->t - p->t_prev) / dlambda };
	const cpp_dec_float_25 rdot{ (p->r - p->r_prev) / dlambda };
	const cpp_dec_float_25 phidot{ (p->phi - p->phi_prev) / dlambda };
	const cpp_dec_float_25 thetadot{ (p->theta - p->theta_prev) / dlambda };

	const cpp_dec_float_25 tsdot{ (ps->t - ps->t_prev) / dlambda };
	const cpp_dec_float_25 rsdot{ (ps->r - ps->r_prev) / dlambda };
	const cpp_dec_float_25 phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const cpp_dec_float_25 thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const cpp_dec_float_25 tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};

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

	const cpp_dec_float_25 temp20{ BH->a8*ps->epsilon };
	const cpp_dec_float_25 temp21{ ps->sigma6 - 0.3333333333333333*BH->a6*ps->epsilon };
	const cpp_dec_float_25 temp22{ BH->a4*ps->epsilon*BH->l4 };
	const cpp_dec_float_25 temp23{ nmag7*p->r3*rdot };
	const cpp_dec_float_25 temp24{ BH->a4*ps->epsilon*BH->l2 };
	const cpp_dec_float_25 temp25{ BH->a4*ps->epsilon };
	const cpp_dec_float_25 temp26{ nmag*p->r*rdot*ps->sigma12 + 0.3333333333333333*BH->a6*ps->epsilon };
	const cpp_dec_float_25 temp27{ nmag*p->r3*rdot };
	const cpp_dec_float_25 temp28{ nmag*p->r5*rdot };
	const cpp_dec_float_25 temp29{ p->mass*nmag13*phidot*p->r*rdot*delta_p + 0.25*BH->a5 };
	const cpp_dec_float_25 temp210{ BH->l4*p->mass*nmag13*phidot*p->r*rdot };
	const cpp_dec_float_25 temp211{ p->mass*nmag13*phidot*p->r3*rdot };
	const cpp_dec_float_25 temp212{ temp211*delta_p + 0.25 };
	const cpp_dec_float_25 temp213{ BH->a3*p->mass*nmag13*phidot };
	const cpp_dec_float_25 temp214{ BH->a*BH->l2*p->mass*nmag13*phidot };
	const cpp_dec_float_25 temp215{ BH->a*p->mass*nmag13 };
	const cpp_dec_float_25 temp216{ BH->a6*ps->epsilon*nmag7*p->r*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp217{ p->r*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp218{ BH->a2*ps->epsilon*BH->l4*nmag7*temp217 };
	const cpp_dec_float_25 temp219{ temp23*ps->sigma6 };
	const cpp_dec_float_25 temp220{ BH->a2*ps->epsilon*BH->l2*temp219 };
	const cpp_dec_float_25 temp221{ BH->a2*ps->epsilon*nmag7*p->r5*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp222{ nmag*p->r*rdot*ps->sigma12*delta_p - 0.3333333333333333 };
	const cpp_dec_float_25 temp223{ BH->a2*ps->epsilon*BH->l4 };
	const cpp_dec_float_25 temp224{ delta_p - 0.3333333333333333*BH->a2 };
	const cpp_dec_float_25 temp225{ delta_p - 0.16666666666666666*BH->a2 };
	const cpp_dec_float_25 temp226{ delta_p - 0.16666666666666666*BH->a6 };
	const cpp_dec_float_25 temp227{ BH->a4*BH->l2*p->mass*nmag13*phidot*p->r*rdot*chi_p };
	const cpp_dec_float_25 temp228{ BH->a4*temp211*chi_p };
	const cpp_dec_float_25 temp229{ p->mass*nmag13*phidot*p->r5*rdot*chi_p };
	const cpp_dec_float_25 temp230{ ps->epsilon*nmag7*temp217*chi_p };
	const cpp_dec_float_25 temp231{ BH->a3*ps->epsilon*BH->l2*nmag7*temp217*chi_p };
	const cpp_dec_float_25 temp232{ chi_p*delta_p - 0.3333333333333333 };
	const cpp_dec_float_25 temp233{ p->r*rdot*ps->sigma12*temp232*BH->a3*ps->epsilon };
	const cpp_dec_float_25 temp234{ p->r*rdot*chi_p2*delta_p + 0.08333333333333333*BH->a3 };
	const cpp_dec_float_25 temp235{ BH->a4*p->mass*nmag13 };
	const cpp_dec_float_25 temp236{ chi_p*delta_p2 - 0.16666666666666666*BH->a2 };
	const cpp_dec_float_25 temp237{ chi_p*delta_p2 - 0.08333333333333333 };
	const cpp_dec_float_25 temp238{ BH->l2*temp211 };
	const cpp_dec_float_25 temp239{ delta_p2 - 0.3333333333333333*BH->a*ps->epsilon };
	const cpp_dec_float_25 temp240{ chi_p*delta_p2 + 0.3333333333333333 };
	const cpp_dec_float_25 temp241{ ps->sigma12*temp240*BH->a*ps->epsilon };
	const cpp_dec_float_25 temp242{ ps->sigma12*chi_p };
	const cpp_dec_float_25 temp243{ temp213*p->r*rdot*chi_p2 };
	const cpp_dec_float_25 temp244{ temp214*p->r*rdot*chi_p2 };
	const cpp_dec_float_25 temp245{ BH->a*temp211*chi_p2 };
	const cpp_dec_float_25 temp246{ BH->a2*ps->epsilon*nmag7*temp217*chi_p2 };
	const cpp_dec_float_25 temp247{ ps->epsilon*nmag*p->r*rdot*ps->sigma12*chi_p2 };
	const cpp_dec_float_25 temp248{ nmag13*phidot*p->r*rdot*chi_p3 };
	const cpp_dec_float_25 temp249{ ps->epsilon*nmag7*temp217*chi_p2 };
	const cpp_dec_float_25 temp250{ BH->a6*ps->epsilon*BH->mass };
	const cpp_dec_float_25 temp251{ sigma_p - 0.3333333333333333*temp24 };
	const cpp_dec_float_25 temp252{ sigma_p - 0.16666666666666666*temp223 };
	const cpp_dec_float_25 temp253{ sigma_p + 0.3333333333333333*temp24 };
	const cpp_dec_float_25 temp254{ sigma_p - 0.3333333333333333*temp25 };
	const cpp_dec_float_25 temp255{ sigma_p - 0.3333333333333333*BH->a2*ps->epsilon*BH->l2 };
	const cpp_dec_float_25 temp256{ temp25*temp219 };
	const cpp_dec_float_25 temp257{ sigma_p - 0.16666666666666666*BH->a2*ps->epsilon };
	const cpp_dec_float_25 temp258{ nmag*rdot*ps->sigma12 };
	const cpp_dec_float_25 temp259{ BH->mass*temp258 };
	const cpp_dec_float_25 temp260{ BH->a6*ps->epsilon };
	const cpp_dec_float_25 temp261{ nmag*p->r*rdot*ps->sigma12 };
	const cpp_dec_float_25 temp262{ nmag*p->r2*rdot*ps->sigma12 };
	const cpp_dec_float_25 temp263{ BH->a2*ps->epsilon*BH->l2 };
	const cpp_dec_float_25 temp264{ temp27*ps->sigma12 };
	const cpp_dec_float_25 temp265{ BH->a2*ps->epsilon*BH->mass };
	const cpp_dec_float_25 temp266{ temp28*ps->sigma12 };
	const cpp_dec_float_25 temp267{ BH->a6*p->mass*BH->mass*nmag13 };
	const cpp_dec_float_25 temp268{ BH->a4*BH->l2*p->mass*BH->mass*nmag13 };
	const cpp_dec_float_25 temp269{ BH->a2*BH->l4*p->mass*BH->mass*nmag13 };
	const cpp_dec_float_25 temp270{ BH->a6*p->mass*nmag13 };
	const cpp_dec_float_25 temp271{ BH->l2*p->mass*nmag13*p->r*rdot*tdot };
	const cpp_dec_float_25 temp272{ p->mass*nmag13*p->r*rdot*tdot };
	const cpp_dec_float_25 temp273{ BH->a4*p->mass*BH->mass*nmag13 };
	const cpp_dec_float_25 temp274{ BH->a2*BH->l2*p->mass*BH->mass*nmag13 };
	const cpp_dec_float_25 temp275{ p->r3*rdot*tdot };
	const cpp_dec_float_25 temp276{ BH->l2*p->mass*nmag13 };
	const cpp_dec_float_25 temp277{ rdot*tdot*sigma_p - 0.08333333333333333 };
	const cpp_dec_float_25 temp278{ phidot*rdot*chi_p };
	const cpp_dec_float_25 temp279{ phidot*p->r*rdot };
	const cpp_dec_float_25 temp280{ chi_p*sigma_p - 0.16666666666666666 };
	const cpp_dec_float_25 temp281{ phidot*p->r2*rdot };
	const cpp_dec_float_25 temp282{ chi_p*sigma_p - 0.08333333333333333 };
	const cpp_dec_float_25 temp283{ chi_p*sigma_p + 0.08333333333333333 };
	const cpp_dec_float_25 temp284{ BH->a5*ps->epsilon*BH->mass };
	const cpp_dec_float_25 temp285{ sigma_p + 0.3333333333333333*BH->a3*ps->epsilon*BH->l2 };
	const cpp_dec_float_25 temp286{ BH->a5*temp230 };
	const cpp_dec_float_25 temp287{ sigma_p + 0.3333333333333333*BH->a3*ps->epsilon };
	const cpp_dec_float_25 temp288{ sigma_p - 0.3333333333333333*BH->a3*ps->epsilon };
	const cpp_dec_float_25 temp289{ nmag*rdot*temp242 };
	const cpp_dec_float_25 temp290{ BH->a5*ps->epsilon*nmag };
	const cpp_dec_float_25 temp291{ p->r2*rdot*temp242 };
	const cpp_dec_float_25 temp292{ sigma_p - 0.08333333333333333*BH->a5*p->mass };
	const cpp_dec_float_25 temp293{ BH->l2*p->mass*BH->mass*nmag13*rdot*tdot };
	const cpp_dec_float_25 temp294{ temp282*BH->a3 };
	const cpp_dec_float_25 temp295{ nmag13*p->r2*rdot*tdot };
	const cpp_dec_float_25 temp296{ temp283*BH->a5 };
	const cpp_dec_float_25 temp297{ nmag13*phidot*rdot*chi_p2 };
	const cpp_dec_float_25 temp298{ BH->l2*p->mass*BH->mass*temp297 };
	const cpp_dec_float_25 temp299{ sigma_p - 0.08333333333333333*BH->a3 };
	const cpp_dec_float_25 temp2100{ BH->a3*p->mass*BH->mass };
	const cpp_dec_float_25 temp2101{ temp211*chi_p2 };
	const cpp_dec_float_25 temp2102{ BH->mass*nmag13*rdot*tdot };
	const cpp_dec_float_25 temp2103{ chi_p*delta_p*sigma_p + 0.08333333333333333 };
	const cpp_dec_float_25 temp2104{ chi_p*delta_p*sigma_p - 0.08333333333333333*BH->a };
	const cpp_dec_float_25 temp2105{ sigma_p + 0.08333333333333333*BH->a };
	const cpp_dec_float_25 temp2106{ delta_p*sigma_p - 0.08333333333333333 };
	const cpp_dec_float_25 temp2107{ delta_p*sigma_p - 0.3333333333333333 };
	const cpp_dec_float_25 temp2108{ delta_p*sigma_p + 0.3333333333333333 };
	const cpp_dec_float_25 temp2109{ temp259*chi_p2 };
	const cpp_dec_float_25 temp2110{ BH->a2*temp247 };
	const cpp_dec_float_25 temp2111{ BH->a2*p->mass };
	const cpp_dec_float_25 temp2112{ chi_p2*temp2106 };
	const cpp_dec_float_25 temp2113{ delta_p*sigma_p + 0.08333333333333333*temp2111 };
	const cpp_dec_float_25 temp2114{ ps->sigma6*chi_p2*delta_p2 };
	const cpp_dec_float_25 temp2115{ ps->epsilon*temp2109 };
	const cpp_dec_float_25 temp2116{ nmag5*rdot*ps->sigma6*sigma_p*(nz*costheta + nx*cosphi*sintheta + ny*sinphi*sintheta) - BH->a6*ps->epsilon };
	const cpp_dec_float_25 temp2117{ nmag5*p->r2*rdot*ps->sigma6*sigma_p };
	const cpp_dec_float_25 temp2118{ ps->epsilon*BH->l2 };
	const cpp_dec_float_25 temp2119{  nx*cosphi*sintheta  };
	const cpp_dec_float_25 temp2120{ p->r4*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp2121{ ps->epsilon*rdot*ps->sigma12 };
	const cpp_dec_float_25 temp2122{ temp260*BH->l2 };
	const cpp_dec_float_25 temp2123{ costheta +temp2119+ ny*sinphi };
	const cpp_dec_float_25 temp2124{ p->r2*rdot*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2125{ p->r4*rdot*ps->sigma12 };
	const cpp_dec_float_25 temp2126{ nmag5*rdot*ps->sigma6*delta_p*sigma_p*(nz*temp2123*sintheta) + BH->a4 };
	const cpp_dec_float_25 temp2127{ ps->epsilon*nmag5*p->r2*rdot*ps->sigma6 };
	const cpp_dec_float_25 temp2128{ temp2118*nmag5*p->r2 };
	const cpp_dec_float_25 temp2129{ delta_p*sigma_p*(nz*temp2123*sintheta) + 0.5 };
	const cpp_dec_float_25 temp2130{ delta_p*sigma_p };
	const cpp_dec_float_25 temp2131{ rdot*ps->sigma12*temp2130 };
	const cpp_dec_float_25 temp2132{ ps->epsilon*BH->l4 };
	const cpp_dec_float_25 temp2133{ p->r2*temp2131 };
	const cpp_dec_float_25 temp2134{ ps->epsilon*temp2125 };
	const cpp_dec_float_25 temp2135{ ps->epsilon*nmag5*rdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2136{ temp2118*nmag5*rdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2137{ temp2127*chi_p };
	const cpp_dec_float_25 temp2138{ BH->a5*ps->epsilon };
	const cpp_dec_float_25 temp2139{ BH->a3*temp2118 };
	const cpp_dec_float_25 temp2140{ BH->a3*ps->epsilon };
	const cpp_dec_float_25 temp2141{ delta_p2*sigma_p*(nz*temp2123*sintheta) - BH->a };
	const cpp_dec_float_25 temp2142{ rdot*temp242*delta_p2*sigma_p };
	const cpp_dec_float_25 temp2143{ BH->a*temp2118 };
	const cpp_dec_float_25 temp2144{ BH->a*ps->epsilon };
	const cpp_dec_float_25 temp2145{ BH->a2*ps->epsilon*nmag5 };
	const cpp_dec_float_25 temp2146{ chi_p2*delta_p2*sigma_p };
	const cpp_dec_float_25 temp2147{ ps->epsilon*nmag5 };
	const cpp_dec_float_25 temp2148{ chi_p2*delta_p3*sigma_p };
	const cpp_dec_float_25 temp2149{ rdot*ps->sigma12 };
	const cpp_dec_float_25 temp2150{ nmag5*phidot*ps->sigma6*sigma_p*(ny*p->r*cosphi*sintheta - nx*p->r*sinphi*sintheta) + BH->a6 };
	const cpp_dec_float_25 temp2151{ temp2147*phidot*p->r2*ps->sigma6 };
	const cpp_dec_float_25 temp2152{ temp2118*nmag5*phidot*p->r2*ps->sigma6 };
	const cpp_dec_float_25 temp2153{ phidot*p->r4*ps->sigma6 };
	const cpp_dec_float_25 temp2154{ cosphi*sintheta - nx*p->r*sinphi };
	const cpp_dec_float_25 temp2155{ phidot*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2156{ phidot*p->r2*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2157{ ps->epsilon*phidot*p->r4*ps->sigma12 };
	const cpp_dec_float_25 temp2158{ temp2147*ps->sigma6*tdot };
	const cpp_dec_float_25 temp2159{ temp2118*nmag5*ps->sigma6*tdot };
	const cpp_dec_float_25 temp2160{ nmag5*p->r2*ps->sigma6*tdot };
	const cpp_dec_float_25 temp2161{ ps->epsilon*ps->sigma12*tdot };
	const cpp_dec_float_25 temp2162{ temp2118*ps->sigma12*tdot };
	const cpp_dec_float_25 temp2163{ ps->epsilon*p->r2*ps->sigma12*tdot };
	const cpp_dec_float_25 temp2164{ nmag5*phidot*ps->sigma6*temp2130*(ny*p->r*temp2154*sintheta) - BH->a4 };
	const cpp_dec_float_25 temp2165{ temp2130*(ny*p->r*temp2154*sintheta) - 0.5 };
	const cpp_dec_float_25 temp2166{ phidot*ps->sigma12*temp2130 };
	const cpp_dec_float_25 temp2167{ phidot*p->r2*ps->sigma12*temp2130 };
	const cpp_dec_float_25 temp2168{ temp2130*(ny*p->r*temp2154*sintheta) + BH->a3 };
	const cpp_dec_float_25 temp2169{ temp2138*ps->sigma12 };
	const cpp_dec_float_25 temp2170{ temp2139*ps->sigma12 };
	const cpp_dec_float_25 temp2171{ temp2140*p->r2*ps->sigma12 };
	const cpp_dec_float_25 temp2172{ temp2147*phidot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2173{ temp2118*nmag5*phidot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2174{ temp2151*chi_p };
	const cpp_dec_float_25 temp2175{ phidot*temp242*temp2130 };
	const cpp_dec_float_25 temp2176{ temp2140*phidot };
	const cpp_dec_float_25 temp2177{ nmag5*ps->sigma6*tdot };
	const cpp_dec_float_25 temp2178{ temp2161*chi_p };
	const cpp_dec_float_25 temp2179{ temp2177*delta_p2*sigma_p*(ny*p->r*temp2154*sintheta) - 0.5 };
	const cpp_dec_float_25 temp2180{ delta_p2*sigma_p };
	const cpp_dec_float_25 temp2181{ temp2180*(ny*p->r*temp2154*sintheta) + BH->a };
	const cpp_dec_float_25 temp2182{ temp242*temp2180 };
	const cpp_dec_float_25 temp2183{ p->r2*temp2182 };
	const cpp_dec_float_25 temp2184{ temp2158*chi_p };
	const cpp_dec_float_25 temp2185{ temp2114*sigma_p };
	const cpp_dec_float_25 temp2186{ phidot*ps->sigma12 };
	const cpp_dec_float_25 temp2187{ delta_p3*sigma_p };
	const cpp_dec_float_25 temp2188{ phidot*ps->sigma6 };
	const cpp_dec_float_25 temp2189{ nmag5*ps->sigma6*thetasdot*sigma_p*(-(nx*ps->r*cosphis*costhetas) - ny*ps->r*costhetas*sinphis + nz*ps->r*sinthetas) - BH->a6 };
	const cpp_dec_float_25 temp2190{ temp2147*p->r2*ps->sigma6*thetasdot };
	const cpp_dec_float_25 temp2191{  ny*ps->r*costhetas*sinphis  };
	const cpp_dec_float_25 temp2192{ thetasdot*sigma_p*(-(nx*ps->r*cosphis*costhetas) -temp2191+ nz*ps->r*sinthetas) - 0.5 };
	const cpp_dec_float_25 temp2193{ (nx*ps->r*cosphis*costhetas) -temp2191 };
	const cpp_dec_float_25 temp2194{ temp2132*ps->sigma12*thetasdot };
	const cpp_dec_float_25 temp2195{ p->r2*ps->sigma12*thetasdot*sigma_p };
	const cpp_dec_float_25 temp2196{ ps->epsilon*p->r4*ps->sigma12*thetasdot };
	const cpp_dec_float_25 temp2197{ nmag5*ps->sigma6*thetasdot*temp2130*(-temp2193+ nz*ps->r*sinthetas) + BH->a4 };
	const cpp_dec_float_25 temp2198{ temp2128*ps->sigma6 };
	const cpp_dec_float_25 temp2199{ temp2145*p->r4 };
	const cpp_dec_float_25 temp2200{ thetasdot*temp2130 };
	const cpp_dec_float_25 temp2201{ ps->epsilon*ps->sigma12 };
	const cpp_dec_float_25 temp2202{ ps->sigma12*temp2200 };
	const cpp_dec_float_25 temp2203{ p->r2*temp2202 };
	const cpp_dec_float_25 temp2204{ temp2147*ps->sigma6*thetasdot*chi_p };
	const cpp_dec_float_25 temp2205{ temp2118*nmag5*ps->sigma6*thetasdot*chi_p };
	const cpp_dec_float_25 temp2206{ temp2190*chi_p };
	const cpp_dec_float_25 temp2207{ thetasdot*chi_p*temp2130 };
	const cpp_dec_float_25 temp2208{ temp2180*(-temp2193+ nz*ps->r*sinthetas) - BH->a };
	const cpp_dec_float_25 temp2209{ ps->sigma12*thetasdot*chi_p*temp2180 };
	const cpp_dec_float_25 temp2210{ temp2144*p->r2 };
	const cpp_dec_float_25 temp2211{ ps->sigma6*thetasdot };
	const cpp_dec_float_25 temp2212{ temp2201*thetasdot };
	const cpp_dec_float_25 temp2213{ ps->sigma12*thetasdot };
	const cpp_dec_float_25 temp2214{ nmag5*rsdot*ps->sigma6*sigma_p*(-(nz*costhetas) - nx*cosphis*sinthetas - ny*sinphis*sinthetas) - BH->a6 };
	const cpp_dec_float_25 temp2215{ temp2147*p->r2*rsdot*ps->sigma6 };
	const cpp_dec_float_25 temp2216{  nx*cosphis*sinthetas  };
	const cpp_dec_float_25 temp2217{ temp2128*rsdot*ps->sigma6 };
	const cpp_dec_float_25 temp2218{ temp25*nmag5*p->r4 };
	const cpp_dec_float_25 temp2219{ (nz*costhetas) -temp2216 };
	const cpp_dec_float_25 temp2220{ rsdot*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2221{ p->r2*temp2220 };
	const cpp_dec_float_25 temp2222{ ps->epsilon*p->r4 };
	const cpp_dec_float_25 temp2223{ nmag5*rsdot*ps->sigma6*temp2130*(-temp2219- ny*sinphis*sinthetas) + BH->a4 };
	const cpp_dec_float_25 temp2224{ temp2130*(-temp2219- ny*sinphis*sinthetas) + 0.5 };
	const cpp_dec_float_25 temp2225{ rsdot*ps->sigma12*temp2130 };
	const cpp_dec_float_25 temp2226{ p->r2*temp2225 };
	const cpp_dec_float_25 temp2227{ temp2147*rsdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2228{ temp2118*nmag5*rsdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2229{ temp2215*chi_p };
	const cpp_dec_float_25 temp2230{ rsdot*temp242*temp2130 };
	const cpp_dec_float_25 temp2231{ temp2180*(-temp2219- ny*sinphis*sinthetas) - BH->a };
	const cpp_dec_float_25 temp2232{ rsdot*temp2182 };
	const cpp_dec_float_25 temp2233{ rsdot*ps->sigma12 };
	const cpp_dec_float_25 temp2234{ rsdot*ps->sigma6 };
	const cpp_dec_float_25 temp2235{ nmag5*phisdot*ps->sigma6*sigma_p*(-(ny*ps->r*cosphis*sinthetas) + nx*ps->r*sinphis*sinthetas) - BH->a6 };
	const cpp_dec_float_25 temp2236{ temp2147*phisdot*p->r2*ps->sigma6 };
	const cpp_dec_float_25 temp2237{ temp2118*nmag5*phisdot*p->r2*ps->sigma6 };
	const cpp_dec_float_25 temp2238{ phisdot*p->r4*ps->sigma6 };
	const cpp_dec_float_25 temp2239{ ps->epsilon*phisdot*ps->sigma12 };
	const cpp_dec_float_25 temp2240{ phisdot*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2241{ (ny*ps->r*cosphis*sinthetas)  };
	const cpp_dec_float_25 temp2242{ p->r2*ps->sigma12*sigma_p };
	const cpp_dec_float_25 temp2243{ ps->epsilon*phisdot*p->r4*ps->sigma12 };
	const cpp_dec_float_25 temp2244{ nmag5*phisdot*ps->sigma6*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) + BH->a4 };
	const cpp_dec_float_25 temp2245{ temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) + 0.5 };
	const cpp_dec_float_25 temp2246{ temp24*phisdot };
	const cpp_dec_float_25 temp2247{ phisdot*ps->sigma12 };
	const cpp_dec_float_25 temp2248{ phisdot*p->r2*ps->sigma12*temp2130 };
	const cpp_dec_float_25 temp2249{ temp2147*phisdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2250{ temp2118*nmag5*phisdot*ps->sigma6*chi_p };
	const cpp_dec_float_25 temp2251{ temp2236*chi_p };
	const cpp_dec_float_25 temp2252{ phisdot*temp242*temp2130 };
	const cpp_dec_float_25 temp2253{ temp2140*phisdot };
	const cpp_dec_float_25 temp2254{ temp2180*(-temp2241+ nx*ps->r*sinphis*sinthetas) - BH->a };
	const cpp_dec_float_25 temp2255{ phisdot*ps->sigma6 };
	const cpp_dec_float_25 temp2256{ delta_p + 2.*BH->a2 };
	const cpp_dec_float_25 temp2257{ chi_p*delta_p - 2.*BH->a };


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((24.*(0. - 0.16666666666666666*temp20*nmag7*p->r*rdot*temp21*BH->l2*nmag7*p->r*rdot*ps->sigma6 - 0.16666666666666666*temp22*nmag7*p->r*rdot*temp21*temp23*ps->sigma6 - 0.3333333333333333*temp24*temp23*ps->sigma6 - 0.16666666666666666*temp25*nmag7*p->r5*rdot*ps->sigma6 + 0.16666666666666666*temp20*temp26*BH->l2*nmag*p->r*rdot*ps->sigma12 + 0.16666666666666666*temp22*temp26*temp27*ps->sigma12 + 0.3333333333333333*temp24*temp27*ps->sigma12 + 0.16666666666666666*temp25*temp28*ps->sigma12 + 0.08333333333333333*BH->a7*temp29*temp276*temp279*delta_p + 0.25*BH->a3*temp210*delta_p + 0.08333333333333333*BH->a*BH->l6*temp29*temp211*delta_p + 0.5*BH->a3*BH->l2*temp212*BH->a*BH->l4*temp212*temp213*p->r5*rdot*delta_p + 0.25*temp214*p->r5*rdot*delta_p + 0.08333333333333333*temp215*phidot*p->r7*rdot*delta_p + 0.16666666666666666*temp216*delta_p + 0.3333333333333333*temp24*nmag7*temp217*delta_p + 0.16666666666666666*temp218*delta_p + 0.3333333333333333*temp256*delta_p + 0.3333333333333333*temp220*delta_p + 0.16666666666666666*temp221*temp226*ps->epsilon*temp222*temp24*temp261*delta_p - 0.16666666666666666*temp223*temp222*temp25*temp264*temp224*temp2118*temp264*temp225*ps->epsilon*temp266*temp226*p->mass*nmag13*temp279*temp232*temp227*temp225*temp210*temp232*temp228*temp224*temp238*chi_p*temp225*temp229*delta_p + 0.3333333333333333*temp286*delta_p + 0.3333333333333333*temp231*delta_p + 0.3333333333333333*temp2140*temp219*temp232*temp290*temp233*BH->l2*nmag*temp233*temp27*temp242*delta_p + 0.08333333333333333*BH->a5*p->mass*nmag13*phidot*temp234*temp276*phidot*temp234*temp2101*delta_p - 0.08333333333333333*temp235*temp279*temp236*temp276*temp279*temp237*temp210*temp236*temp211*chi_p*delta_p2 - 0.16666666666666666*temp238*temp237*temp229*delta_p2 - 0.3333333333333333*BH->a3*temp230*temp239*BH->l2*nmag7*temp217*chi_p*temp239*temp219*temp240*temp2140*nmag*p->r*rdot*temp241*BH->l2*nmag*p->r*rdot*temp241*temp27*temp242*delta_p2 + 0.16666666666666666*temp243*delta_p2 + 0.16666666666666666*temp244*delta_p2 + 0.16666666666666666*temp245*delta_p2 - 0.16666666666666666*temp246*delta_p2 + 0.16666666666666666*temp2110*delta_p2 - 0.08333333333333333*temp2111*temp248*delta_p2 + 0.16666666666666666*temp249*delta_p3 - 0.16666666666666666*temp247*delta_p3 - 0.16666666666666666*temp250*nmag7*rdot*ps->sigma6*temp251*BH->mass*nmag7*rdot*ps->sigma6*temp252*BH->mass*nmag7*rdot*ps->sigma6*sigma_p + 0.16666666666666666*temp216*temp253*nmag7*temp217*sigma_p + 0.16666666666666666*temp218*temp254*BH->mass*nmag7*p->r2*rdot*ps->sigma6*temp255*BH->mass*nmag7*p->r2*rdot*ps->sigma6*sigma_p + 0.3333333333333333*temp256*sigma_p + 0.3333333333333333*temp220*temp257*BH->mass*nmag7*temp2120*sigma_p + 0.16666666666666666*temp221*sigma_p + 0.16666666666666666*temp250*temp258*temp253*temp259*sigma_p + 0.16666666666666666*temp223*temp259*sigma_p - 0.16666666666666666*temp260*temp261*temp251*temp261*temp252*temp261*sigma_p + 0.3333333333333333*temp25*BH->mass*temp262*sigma_p + 0.3333333333333333*temp263*BH->mass*temp262*temp254*temp264*temp255*temp264*sigma_p + 0.16666666666666666*temp265*nmag*temp2125*temp257*temp266*sigma_p + 0.08333333333333333*temp267*rdot*tdot*sigma_p + 0.16666666666666666*temp268*rdot*tdot*sigma_p + 0.08333333333333333*temp269*temp277*temp270*p->r*rdot*tdot*sigma_p - 0.16666666666666666*BH->a4*temp271*sigma_p - 0.08333333333333333*BH->a2*BH->l4*temp272*sigma_p + 0.16666666666666666*temp273*p->r2*rdot*tdot*sigma_p + 0.16666666666666666*temp274*p->r2*rdot*tdot*sigma_p - 0.16666666666666666*temp235*temp275*sigma_p - 0.16666666666666666*BH->a2*temp276*temp275*sigma_p + 0.08333333333333333*temp2111*BH->mass*nmag13*p->r4*temp277*temp2111*nmag13*p->r5*temp277*temp267*temp278*sigma_p - 0.16666666666666666*temp268*temp278*sigma_p - 0.08333333333333333*temp269*temp278*sigma_p + 0.08333333333333333*temp270*temp279*chi_p*sigma_p + 0.16666666666666666*temp227*sigma_p + 0.08333333333333333*BH->a2*temp210*temp280*temp273*temp281*temp280*temp274*temp281*chi_p*sigma_p + 0.16666666666666666*temp228*sigma_p + 0.16666666666666666*BH->a2*temp238*temp282*temp2111*BH->mass*nmag13*phidot*p->r4*rdot*temp283*BH->a2*temp229*sigma_p + 0.3333333333333333*temp284*nmag7*rdot*ps->sigma6*chi_p*temp285*BH->mass*nmag7*rdot*ps->sigma6*chi_p*sigma_p - 0.3333333333333333*temp286*sigma_p - 0.3333333333333333*temp231*temp287*BH->mass*nmag7*p->r2*rdot*ps->sigma6*chi_p*temp288*temp219*chi_p*sigma_p - 0.3333333333333333*temp284*temp289*temp288*BH->l2*BH->mass*temp289*sigma_p + 0.3333333333333333*temp290*p->r*rdot*temp242*temp285*nmag*p->r*rdot*temp242*temp288*BH->mass*nmag*temp291*temp287*temp27*temp242*temp292*temp2102*temp294*temp293*temp296*temp272*temp283*BH->a3*temp271*temp294*p->mass*BH->mass*temp295*temp283*BH->a3*p->mass*nmag13*temp275*temp296*p->mass*BH->mass*temp297*sigma_p + 0.08333333333333333*BH->a3*temp298*temp292*nmag13*temp279*chi_p2*temp299*temp276*temp279*chi_p2*sigma_p + 0.08333333333333333*temp2100*nmag13*temp281*chi_p2*temp299*temp2101*temp299*p->mass*temp2102*temp2104*temp293*temp2103*BH->a3*temp272*temp2103*BH->a*temp271*temp2104*p->mass*BH->mass*temp295*temp2103*temp215*temp275*temp2103*temp2100*temp297*delta_p*temp2105*temp298*temp2106*temp243*temp2106*temp244*delta_p*temp2105*p->mass*BH->mass*nmag13*temp281*temp2112*temp245*temp2107*temp265*nmag7*rdot*ps->sigma6*chi_p2*temp2108*temp246*temp2108*BH->a2*temp2115*temp2107*temp2110*temp2113*temp2102*temp2112*BH->a2*temp272*temp2112*temp2111*BH->mass*nmag13*phidot*rdot*chi_p3*temp2113*temp248*delta_p*sigma_p + 0.16666666666666666*ps->epsilon*BH->mass*nmag7*rdot*temp2114*sigma_p - 0.16666666666666666*temp249*delta_p2*sigma_p - 0.16666666666666666*temp2115*delta_p2*sigma_p + 0.16666666666666666*temp247*delta_p2*sigma_p - 0.5*temp20*temp2116*BH->l2*nmag5*rdot*ps->sigma6*sigma_p*(nz*temp2123*sintheta) - 0.5*temp22*temp2116*temp2117*(nz*temp2123*sintheta) - BH->a4*temp2118*temp2117*(nz*temp2123*sintheta) - 0.5*temp25*nmag5*temp2120*sigma_p*(nz*temp2123*sintheta) + (BH->a8*temp2121*sigma_p*(nz*temp2123*sintheta))/nmag + (2.*temp2122*temp2149*sigma_p*(nz*temp2123*sintheta))/nmag + (BH->a4*temp2132*temp2149*sigma_p*(nz*temp2123*sintheta))/nmag + (2.*temp260*temp2124*(nz*temp2123*sintheta))/nmag + (2.*temp24*temp2124*(nz*temp2123*sintheta))/nmag + (BH->a4*temp2134*sigma_p*(nz*temp2123*sintheta))/nmag + 0.5*temp260*temp2126*temp2118*nmag5*rdot*ps->sigma6*temp2129*temp223*temp2126*temp2127*temp2130*(nz*temp2123*sintheta) + BH->a2*temp2128*rdot*ps->sigma6*temp2129*temp2145*temp2120*temp2130*(nz*temp2123*sintheta) - (BH->a6*temp2121*temp2130*(nz*temp2123*sintheta))/nmag - (2.*temp24*temp2131*(nz*temp2123*sintheta))/nmag - (BH->a2*temp2132*temp2131*(nz*temp2123*sintheta))/nmag - (2.*temp25*temp2133*(nz*temp2123*sintheta))/nmag - (2.*temp263*temp2133*(nz*temp2123*sintheta))/nmag - (BH->a2*temp2134*temp2130*(nz*temp2123*sintheta))/nmag + BH->a5*temp2135*temp2130*(nz*temp2123*sintheta) + BH->a3*temp2136*temp2130*(nz*temp2123*sintheta) + BH->a3*temp2137*temp2130*(nz*temp2123*sintheta) - (2.*temp2138*rdot*temp242*temp2130*(nz*temp2123*sintheta))/nmag - (2.*temp2139*rdot*temp242*temp2130*(nz*temp2123*sintheta))/nmag - (2.*temp2140*temp291*temp2130*(nz*temp2123*sintheta))/nmag - BH->a3*temp2135*temp2141*temp2136*temp2141*temp2137*temp2180*(nz*temp2123*sintheta) + (2.*temp2140*temp2142*(nz*temp2123*sintheta))/nmag + (2.*temp2143*temp2142*(nz*temp2123*sintheta))/nmag + (2.*temp2144*temp291*temp2180*(nz*temp2123*sintheta))/nmag - 0.5*temp2145*rdot*temp2185*(nz*temp2123*sintheta) + (BH->a2*temp2121*temp2146*(nz*temp2123*sintheta))/nmag + 0.5*temp2147*rdot*ps->sigma6*temp2148*(nz*temp2123*sintheta) - (ps->epsilon*temp2149*temp2148*(nz*temp2123*sintheta))/nmag + 0.5*temp20*temp2150*temp2118*nmag5*temp2188*sigma_p*(ny*p->r*temp2154*sintheta) + 0.5*temp22*temp2150*temp2151*sigma_p*(ny*p->r*temp2154*sintheta) + BH->a4*temp2152*sigma_p*(ny*p->r*temp2154*sintheta) + 0.5*temp25*nmag5*temp2153*sigma_p*(ny*p->r*temp2154*sintheta) - (BH->a8*ps->epsilon*temp2155*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp2122*temp2155*(ny*p->r*temp2154*sintheta))/nmag - (BH->a4*temp2132*temp2155*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp260*temp2156*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp24*temp2156*(ny*p->r*temp2154*sintheta))/nmag - (BH->a4*temp2157*sigma_p*(ny*p->r*temp2154*sintheta))/nmag - 0.5*BH->a7*temp2158*sigma_p*(ny*p->r*temp2154*sintheta) - 0.5*BH->a5*temp2159*sigma_p*(ny*p->r*temp2154*sintheta) - 0.5*temp2138*temp2160*sigma_p*(ny*p->r*temp2154*sintheta) + (BH->a7*temp2161*sigma_p*(ny*p->r*temp2154*sintheta))/nmag + (BH->a5*temp2162*sigma_p*(ny*p->r*temp2154*sintheta))/nmag + (BH->a5*temp2163*sigma_p*(ny*p->r*temp2154*sintheta))/nmag - 0.5*temp260*temp2164*temp2118*nmag5*temp2188*temp2165*temp223*temp2164*temp2151*temp2130*(ny*p->r*temp2154*sintheta) - BH->a2*temp2152*temp2165*temp2145*temp2153*temp2130*(ny*p->r*temp2154*sintheta) + (BH->a6*ps->epsilon*temp2166*(ny*p->r*temp2154*sintheta))/nmag + (2.*temp24*temp2166*(ny*p->r*temp2154*sintheta))/nmag + (BH->a2*temp2132*temp2166*(ny*p->r*temp2154*sintheta))/nmag + (2.*temp25*temp2167*(ny*p->r*temp2154*sintheta))/nmag + (2.*temp263*temp2167*(ny*p->r*temp2154*sintheta))/nmag + (BH->a2*temp2157*temp2130*(ny*p->r*temp2154*sintheta))/nmag + BH->a5*temp2158*temp2168*temp2159*temp2168*temp2147*p->r2*ps->sigma6*tdot*temp2130*(ny*p->r*temp2154*sintheta) - (2.*temp2169*tdot*temp2130*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp2170*tdot*temp2130*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp2171*tdot*temp2130*(ny*p->r*temp2154*sintheta))/nmag - BH->a5*temp2172*temp2130*(ny*p->r*temp2154*sintheta) - BH->a3*temp2173*temp2130*(ny*p->r*temp2154*sintheta) - BH->a3*temp2174*temp2130*(ny*p->r*temp2154*sintheta) + (2.*temp2138*temp2175*(ny*p->r*temp2154*sintheta))/nmag + (2.*temp2139*temp2175*(ny*p->r*temp2154*sintheta))/nmag + (2.*temp2176*p->r2*temp242*temp2130*(ny*p->r*temp2154*sintheta))/nmag + 0.5*temp25*temp2177*chi_p*temp2130*(ny*p->r*temp2154*sintheta) - (BH->a4*temp2178*temp2130*(ny*p->r*temp2154*sintheta))/nmag - 0.5*temp2140*temp2179*temp2143*temp2179*temp2144*temp2160*temp2180*(ny*p->r*temp2154*sintheta) + (BH->a3*temp2161*temp2180*(ny*p->r*temp2154*sintheta))/nmag + (BH->a*temp2162*temp2180*(ny*p->r*temp2154*sintheta))/nmag + (BH->a*temp2163*temp2180*(ny*p->r*temp2154*sintheta))/nmag + BH->a3*temp2172*temp2181*temp2173*temp2181*temp2174*temp2180*(ny*p->r*temp2154*sintheta) - (2.*temp2176*temp2182*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp2143*phidot*temp2182*(ny*p->r*temp2154*sintheta))/nmag - (2.*temp2144*phidot*temp2183*(ny*p->r*temp2154*sintheta))/nmag - BH->a2*temp2184*temp2180*(ny*p->r*temp2154*sintheta) + (2.*BH->a2*temp2178*temp2180*(ny*p->r*temp2154*sintheta))/nmag + 0.5*temp2145*phidot*temp2185*(ny*p->r*temp2154*sintheta) - (BH->a2*ps->epsilon*temp2186*temp2146*(ny*p->r*temp2154*sintheta))/nmag + 0.5*temp2184*temp2187*(ny*p->r*temp2154*sintheta) - (ps->epsilon*ps->sigma12*tdot*chi_p*temp2187*(ny*p->r*temp2154*sintheta))/nmag - 0.5*temp2147*temp2188*temp2148*(ny*p->r*temp2154*sintheta) + (ps->epsilon*temp2186*temp2148*(ny*p->r*temp2154*sintheta))/nmag - 0.5*temp20*temp2189*temp2118*nmag5*ps->sigma6*temp2192*temp22*temp2189*temp2190*sigma_p*(-temp2193+ nz*ps->r*sinthetas) - BH->a4*temp2198*temp2192*temp2218*temp2211*sigma_p*(-temp2193+ nz*ps->r*sinthetas) + (BH->a8*temp2212*sigma_p*(-temp2193+ nz*ps->r*sinthetas))/nmag + (2.*temp2122*temp2213*sigma_p*(-temp2193+ nz*ps->r*sinthetas))/nmag + (BH->a4*temp2194*sigma_p*(-temp2193+ nz*ps->r*sinthetas))/nmag + (2.*temp260*temp2195*(-temp2193+ nz*ps->r*sinthetas))/nmag + (2.*temp24*temp2195*(-temp2193+ nz*ps->r*sinthetas))/nmag + (BH->a4*temp2196*sigma_p*(-temp2193+ nz*ps->r*sinthetas))/nmag + 0.5*temp260*temp2197*temp2118*nmag5*ps->sigma6*temp2200*(-temp2193+ nz*ps->r*sinthetas) + 0.5*temp223*temp2197*temp2190*temp2130*(-temp2193+ nz*ps->r*sinthetas) + BH->a2*temp2198*temp2200*(-temp2193+ nz*ps->r*sinthetas) + 0.5*temp2199*ps->sigma6*temp2200*(-temp2193+ nz*ps->r*sinthetas) - (BH->a6*temp2201*temp2200*(-temp2193+ nz*ps->r*sinthetas))/nmag - (2.*temp24*temp2202*(-temp2193+ nz*ps->r*sinthetas))/nmag - (BH->a2*temp2194*temp2130*(-temp2193+ nz*ps->r*sinthetas))/nmag - (2.*temp25*temp2203*(-temp2193+ nz*ps->r*sinthetas))/nmag - (2.*temp263*temp2203*(-temp2193+ nz*ps->r*sinthetas))/nmag - (BH->a2*temp2196*temp2130*(-temp2193+ nz*ps->r*sinthetas))/nmag + BH->a5*temp2204*temp2130*(-temp2193+ nz*ps->r*sinthetas) + BH->a3*temp2205*temp2130*(-temp2193+ nz*ps->r*sinthetas) + BH->a3*temp2206*temp2130*(-temp2193+ nz*ps->r*sinthetas) - (2.*temp2169*temp2207*(-temp2193+ nz*ps->r*sinthetas))/nmag - (2.*temp2170*temp2207*(-temp2193+ nz*ps->r*sinthetas))/nmag - (2.*temp2171*temp2207*(-temp2193+ nz*ps->r*sinthetas))/nmag - BH->a3*temp2204*temp2208*temp2205*temp2208*temp2206*temp2180*(-temp2193+ nz*ps->r*sinthetas) + (2.*temp2140*temp2209*(-temp2193+ nz*ps->r*sinthetas))/nmag + (2.*temp2143*temp2209*(-temp2193+ nz*ps->r*sinthetas))/nmag + (2.*temp2210*temp2209*(-temp2193+ nz*ps->r*sinthetas))/nmag - 0.5*temp2145*temp2211*temp2146*(-temp2193+ nz*ps->r*sinthetas) + (BH->a2*temp2212*temp2146*(-temp2193+ nz*ps->r*sinthetas))/nmag + 0.5*temp2147*temp2211*temp2148*(-temp2193+ nz*ps->r*sinthetas) - (ps->epsilon*temp2213*temp2148*(-temp2193+ nz*ps->r*sinthetas))/nmag - 0.5*temp20*temp2214*temp2118*nmag5*temp2234*sigma_p*(-temp2219- ny*sinphis*sinthetas) - 0.5*temp22*temp2214*temp2215*sigma_p*(-temp2219- ny*sinphis*sinthetas) - BH->a4*temp2217*sigma_p*(-temp2219- ny*sinphis*sinthetas) - 0.5*temp2218*temp2234*sigma_p*(-temp2219- ny*sinphis*sinthetas) + (BH->a8*ps->epsilon*temp2220*(-temp2219- ny*sinphis*sinthetas))/nmag + (2.*temp2122*temp2220*(-temp2219- ny*sinphis*sinthetas))/nmag + (BH->a4*temp2132*temp2220*(-temp2219- ny*sinphis*sinthetas))/nmag + (2.*temp260*temp2221*(-temp2219- ny*sinphis*sinthetas))/nmag + (2.*temp24*temp2221*(-temp2219- ny*sinphis*sinthetas))/nmag + (BH->a4*temp2222*temp2220*(-temp2219- ny*sinphis*sinthetas))/nmag + 0.5*temp260*temp2223*temp2118*nmag5*temp2234*temp2224*temp223*temp2223*temp2215*temp2130*(-temp2219- ny*sinphis*sinthetas) + BH->a2*temp2217*temp2224*temp2199*temp2234*temp2130*(-temp2219- ny*sinphis*sinthetas) - (BH->a6*ps->epsilon*temp2225*(-temp2219- ny*sinphis*sinthetas))/nmag - (2.*temp24*temp2225*(-temp2219- ny*sinphis*sinthetas))/nmag - (BH->a2*temp2132*temp2225*(-temp2219- ny*sinphis*sinthetas))/nmag - (2.*temp25*temp2226*(-temp2219- ny*sinphis*sinthetas))/nmag - (2.*temp263*temp2226*(-temp2219- ny*sinphis*sinthetas))/nmag - (BH->a2*temp2222*temp2225*(-temp2219- ny*sinphis*sinthetas))/nmag + BH->a5*temp2227*temp2130*(-temp2219- ny*sinphis*sinthetas) + BH->a3*temp2228*temp2130*(-temp2219- ny*sinphis*sinthetas) + BH->a3*temp2229*temp2130*(-temp2219- ny*sinphis*sinthetas) - (2.*temp2138*temp2230*(-temp2219- ny*sinphis*sinthetas))/nmag - (2.*temp2139*temp2230*(-temp2219- ny*sinphis*sinthetas))/nmag - (2.*temp2140*p->r2*temp2230*(-temp2219- ny*sinphis*sinthetas))/nmag - BH->a3*temp2227*temp2231*temp2228*temp2231*temp2229*temp2180*(-temp2219- ny*sinphis*sinthetas) + (2.*temp2140*temp2232*(-temp2219- ny*sinphis*sinthetas))/nmag + (2.*temp2143*temp2232*(-temp2219- ny*sinphis*sinthetas))/nmag + (2.*temp2210*temp2232*(-temp2219- ny*sinphis*sinthetas))/nmag - 0.5*temp2145*rsdot*temp2185*(-temp2219- ny*sinphis*sinthetas) + (BH->a2*ps->epsilon*temp2233*temp2146*(-temp2219- ny*sinphis*sinthetas))/nmag + 0.5*temp2147*temp2234*temp2148*(-temp2219- ny*sinphis*sinthetas) - (ps->epsilon*temp2233*temp2148*(-temp2219- ny*sinphis*sinthetas))/nmag - 0.5*temp20*temp2235*temp2118*nmag5*temp2255*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas) - 0.5*temp22*temp2235*temp2236*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas) - BH->a4*temp2237*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas) - 0.5*temp25*nmag5*temp2238*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas) + (BH->a8*temp2239*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (2.*temp2122*temp2240*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (BH->a4*temp2132*temp2240*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (2.*temp260*phisdot*temp2242*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (2.*temp2246*temp2242*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (BH->a4*temp2243*sigma_p*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + 0.5*temp260*temp2244*temp2118*nmag5*temp2255*temp2245*temp223*temp2244*temp2236*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) + BH->a2*temp2237*temp2245*temp2145*temp2238*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) - (BH->a6*temp2239*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (2.*temp2246*ps->sigma12*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (BH->a2*temp2132*temp2247*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (2.*temp25*temp2248*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (2.*temp263*temp2248*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (BH->a2*temp2243*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + BH->a5*temp2249*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) + BH->a3*temp2250*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) + BH->a3*temp2251*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas) - (2.*temp2138*temp2252*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (2.*temp2139*temp2252*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - (2.*temp2253*p->r2*temp242*temp2130*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - BH->a3*temp2249*temp2254*temp2250*temp2254*temp2251*temp2180*(-temp2241+ nx*ps->r*sinphis*sinthetas) + (2.*temp2253*temp2182*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (2.*temp2143*phisdot*temp2182*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + (2.*temp2144*phisdot*temp2183*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag - 0.5*temp2145*phisdot*temp2185*(-temp2241+ nx*ps->r*sinphis*sinthetas) + (BH->a2*temp2239*temp2146*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag + 0.5*temp2147*temp2255*temp2148*(-temp2241+ nx*ps->r*sinphis*sinthetas) - (ps->epsilon*temp2247*temp2148*(-temp2241+ nx*ps->r*sinphis*sinthetas))/nmag))/(p->mass*nmag13*(BH->a3 + BH->a*BH->l2 + BH->a*p->r2 - chi_p*delta_p)*(BH->a4*temp2256*BH->l2*delta_p + BH->l4*temp2256*p->r2*delta_p + 2.*BH->l2*p->r2*delta_p + p->r4*delta_p - 2.*BH->a3*temp2257*BH->l2*temp2257*p->r2*chi_p*delta_p + BH->a2*chi_p2*delta_p)*sigma_p));
	return phi_new;
}

cpp_dec_float_25 calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 theta_new{};
	const cpp_dec_float_25 sigma_p2{ sigma_p*sigma_p };
	const cpp_dec_float_25 sigma_p3{ sigma_p*sigma_p2 };
	const cpp_dec_float_25 sigma_p4{ sigma_p2*sigma_p2 };
	const cpp_dec_float_25 delta_p2{ delta_p*delta_p };
	const cpp_dec_float_25 delta_p3{ delta_p*delta_p2 };
	const cpp_dec_float_25 delta_p4{ delta_p2*delta_p2 };
	const cpp_dec_float_25 chi_p2{ chi_p*chi_p };
	const cpp_dec_float_25 chi_p3{ chi_p*chi_p2 };
	const cpp_dec_float_25 chi_p4{ chi_p2*chi_p2 };

	const cpp_dec_float_25 tdot{ (p->t - p->t_prev) / dlambda };
	const cpp_dec_float_25 rdot{ (p->r - p->r_prev) / dlambda };
	const cpp_dec_float_25 phidot{ (p->phi - p->phi_prev) / dlambda };
	const cpp_dec_float_25 thetadot{ (p->theta - p->theta_prev) / dlambda };

	const cpp_dec_float_25 tsdot{ (ps->t - ps->t_prev) / dlambda };
	const cpp_dec_float_25 rsdot{ (ps->r - ps->r_prev) / dlambda };
	const cpp_dec_float_25 phisdot{ (ps->phi - ps->phi_prev) / dlambda };
	const cpp_dec_float_25 thetasdot{ (ps->theta - ps->theta_prev) / dlambda };
	const cpp_dec_float_25 tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)};
	const cpp_dec_float_25 thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)};

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



	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*((-0.5*(0. - (96.*ps->epsilon*ps->sigma6*(-0.5*nmag6 + ps->sigma6)*(tdot*(-0.5*BH->a2 + 0.5*delta_p) + phidot*(BH->a3 + BH->a*BH->l2 + BH->a*p->r2 - chi_p*delta_p))*(nx*p->r*cosphi*costheta + ny*p->r*costheta*sinphi - nz*p->r*sintheta))/(p->mass*nmag14*sigma_p)))/sigma_p);
	return theta_new;
}

std::tuple<cpp_dec_float_25, cpp_dec_float_25,
	cpp_dec_float_25, cpp_dec_float_25> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	cpp_dec_float_25 t_new{};
	cpp_dec_float_25 r_new{};
	cpp_dec_float_25 phi_new{};
	cpp_dec_float_25 theta_new{};

	t_new = calc_t_next(BH,
	p, ps,
	sigma_p, delta_p,
	chi_p, dlambda);
	r_new = calc_r_next(BH,
	p, ps,
	sigma_p, delta_p,
	chi_p, dlambda);
	phi_new = calc_phi_next(BH,
	p, ps,
	sigma_p, delta_p,
	chi_p, dlambda);
	theta_new = calc_theta_next(BH,
	p, ps,
	sigma_p, delta_p,
	chi_p, dlambda);
	return {t_new, r_new, phi_new, theta_new};
}
