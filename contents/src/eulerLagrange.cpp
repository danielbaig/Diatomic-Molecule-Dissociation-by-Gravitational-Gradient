// Configured==true
#include "eulerLagrange.h"

const cpp_dec_float_n epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


Precalculated::Precalculated(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda):
	sigma_p2{ sigma_p*sigma_p },
	sigma_p3{ sigma_p*sigma_p2 },
	sigma_p4{ sigma_p2*sigma_p2 },
	delta_p2{ delta_p*delta_p },
	delta_p3{ delta_p*delta_p2 },
	delta_p4{ delta_p2*delta_p2 },
	chi_p2{ chi_p*chi_p },
	chi_p3{ chi_p*chi_p2 },
	chi_p4{ chi_p2*chi_p2 },

	tdot{ (p->t_prev - p->t_prevprev) / (dlambda) },
	rdot{ (p->r_prev - p->r_prevprev) / (dlambda) },
	phidot{ (p->phi_prev - p->phi_prevprev) / (dlambda) },
	thetadot{ (p->theta_prev - p->theta_prevprev) / (dlambda) },

	tsdot{ (ps->t_prev - ps->t_prevprev) / (dlambda) },
	rsdot{ (ps->r_prev - ps->r_prevprev) / (dlambda) },
	phisdot{ (ps->phi_prev - ps->phi_prevprev) / (dlambda) },
	thetasdot{ (ps->theta_prev - ps->theta_prevprev) / (dlambda) },
	tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)},
	rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)},
	phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)},
	thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)},

	nx{p->r*sin(p->theta)*cos(p->phi) - ps->r*sin(ps->theta)*cos(ps->phi)},
	ny{p->r*sin(p->theta)*sin(p->phi) - ps->r*sin(ps->theta)*sin(ps->phi)},
	nz{p->r*cos(p->theta) - ps->r*cos(ps->theta)},
	nmag{sqrt(nx*nx + ny*ny + nz*nz)},
	nmag2{nmag*nmag},
	nmag5{nmag*nmag2*nmag2},
	nmag6{nmag5*nmag},
	nmag7{nmag6*nmag},
	nmag8{nmag6*nmag2},
	nmag12{nmag6*nmag6},
	nmag13{nmag6*nmag7},
	nmag14{nmag7*nmag7},

	ral_sqr{BH->a2 + BH->l2 + p->r2},
	sintheta{sin(p->theta)},
	sintheta2{sintheta*sintheta},
	sintheta3{sintheta2*sintheta},
	sintheta4{sintheta2*sintheta2},
	sintheta5{sintheta2*sintheta3},
	sintheta6{sintheta3*sintheta3},
		sintheta7{sintheta3*sintheta4},

	costheta{cos(p->theta)},
	sin2theta{sin(2.*p->theta)},
	sinphi{sin(p->phi)},
	cosphi{cos(p->phi)},

	sinthetas{sin(ps->theta)},
	costhetas{cos(ps->theta)},
	sinphis{sin(ps->phi)},
	cosphis{cos(ps->phi)},
	sinthetasdot{sin(thetasdot)},
	costhetasdot{cos(thetasdot)},
	sinphisdot{sin(phisdot)},
	cosphisdot{cos(phisdot)},
	sinthetas2{sinthetas*sinthetas},
	costhetas2{costhetas*costhetas},
	sinphis2{sinphis*sinphis},
	cosphis2{cosphis*cosphis},
	rsdot2{rsdot*rsdot},
	phisdot2{phisdot*phisdot}
{
	/*
	Initalise the struct.
	*/
	return;
}


Precalculated::~Precalculated()
{
	/*
	Delete the struct.
	*/
delete[] temp0;
temp0 = nullptr;
delete[] temp1;
temp1 = nullptr;
delete[] temp2;
temp2 = nullptr;
}


cpp_dec_float_n calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of t.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n t_new: New coordinate.
	*/

	cpp_dec_float_n t_new{};
	temp[0] = -1.*BH->a4-2.*BH->a2*p->r2-p->r4+BH->a2*delta_p; // (...)
	temp[1] = BH->charge*p->r; // *...+
	temp[2] = precalculated->phidot*p->r2; // *...+
	temp[3] = -1.*BH->a2+delta_p; // (...)
	temp[4] = BH->a3-BH->a*delta_p; // (...)
	temp[5] = 2.*precalculated->rdot*(p->charge*BH->charge-2.*BH->a*p->mass*precalculated->phidot*p->r+(-2.*BH->mass+2.*p->r)*(BH->a*p->mass*precalculated->phidot-p->mass*precalculated->tdot)); // (...)
	temp[6] = -p->charge*BH->charge; // -...*
	temp[7] = precalculated->tdot*temp[3]+precalculated->phidot*temp[4]; // (...)
	temp[8] = -4.*p->mass*precalculated->phidot; // -...*
	temp[9] = BH->a*precalculated->phidot-precalculated->tdot; // (...)
	temp[10] = p->charge*BH->charge; // +...*
	temp[11] = -4.*BH->a2*(BH->a2+p->r2-delta_p)*(-1.*BH->a2-p->r2+delta_p); // (...)
	temp[12] = -4.*p->r*precalculated->rdot*(-1.*p->charge*temp[1]+BH->a*p->mass*temp[2]+p->mass*precalculated->tdot*temp[3]+p->mass*precalculated->phidot*temp[4]); // (...)
	temp[13] = precalculated->phidot*p->r4; // *...+
	temp[14] = 2.*BH->a*precalculated->phidot-precalculated->tdot; // (...)
	temp[15] = 2.*precalculated->rdot*(+temp[8]*p->r3+BH->a*p->mass*(-2.*BH->mass+2.*p->r)*temp[9]+BH->a*p->mass*p->r*(-4.*BH->a*precalculated->phidot+2.*precalculated->tdot)+temp[10]*chi_p); // (...)
	temp[16] = -2.*temp[0]*(temp[12]/(p->mass*precalculated->sigma_p2)-temp[5]/(p->mass*sigma_p)); // (...)
	temp[17] = BH->a*p->mass; // +...*
	temp[18] = temp[6]*p->r*chi_p; // +...+
	temp[19] = 4.*(BH->a2-delta_p)*temp[0]; // (...)
	temp[20] = -1.*BH->a2-p->r2+delta_p; // (...)
	temp[21] = p->mass*temp[13]+temp[17]*p->r2*temp[14]+temp[18]+temp[17]*temp[7]; // (...)
	temp[22] = BH->a*temp[20]; // *...*
	temp[23] = (4.*p->r*precalculated->rdot*temp[21])/(p->mass*precalculated->sigma_p2)+temp[15]/(p->mass*sigma_p); // (...)
	temp[24] = temp[16]/sigma_p-(2.*temp[22]*temp[23])/sigma_p; // (...)


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-(temp[24]/(temp[11]/precalculated->sigma_p2+temp[19]/precalculated->sigma_p2)));
	return t_new;
}

cpp_dec_float_n calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of r.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n r_new: New coordinate.
	*/

	cpp_dec_float_n r_new{};
	temp[0] = -2.*p->r*(precalculated->rdot*precalculated->rdot); // (...)
	temp[1] = (4.*p->charge*BH->charge*p->r2*(-1.*precalculated->tdot+precalculated->phidot*chi_p))/(p->mass*precalculated->sigma_p2); // +...+
	temp[2] = -1.*BH->a*precalculated->phidot; // -...+
	temp[3] = precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot; // (...)
	temp[4] = (-2.*BH->mass+2.*p->r)*(precalculated->rdot*precalculated->rdot)*sigma_p; // (...)
	temp[5] = -(2.*p->r*pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2))/precalculated->sigma_p2; // -...+
	temp[6] = (4.*precalculated->phidot*p->r*temp[3])/sigma_p; // +...-
	temp[7] = temp[5]+temp[1]; // +...+
	temp[8] = -((-2.*BH->mass+2.*p->r)*pow(+temp[2]+precalculated->tdot,2))/sigma_p; // -...+
	temp[9] = temp[7]+(2.*p->r*pow(+temp[2]+precalculated->tdot,2)*delta_p)/precalculated->sigma_p2; // +...+
	temp[10] = temp[9]+temp[8]; // +...+
	temp[11] = 2.*p->charge*BH->charge*(-1.*precalculated->tdot+precalculated->phidot*chi_p); // (...)
	temp[12] = temp[10]+temp[6]; // +...-
	temp[13] = temp[12]-temp[11]/(p->mass*sigma_p); // +...+


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((0.5*delta_p*(temp[0]/delta_p+temp[13]+temp[4]/precalculated->delta_p2))/sigma_p);
	return r_new;
}

cpp_dec_float_n calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of phi.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n phi_new: New coordinate.
	*/

	cpp_dec_float_n phi_new{};
	temp[0] = -1.*BH->a5*p->charge; // -...*
	temp[1] = -2.*BH->a3*p->charge; // -...*
	temp[2] = -BH->a*p->charge; // -...*
	temp[3] = p->r6*precalculated->rdot; // *...+
	temp[4] = p->r2*precalculated->rdot; // *...*
	temp[5] = BH->charge*p->r4; // *...*
	temp[6] = 2.*BH->a3*p->charge; // +...*
	temp[7] = 2.*BH->a*p->charge; // +...*
	temp[8] = BH->a2*p->mass; // +...*
	temp[9] = p->r5*precalculated->rdot; // *...*
	temp[10] = p->r7*precalculated->rdot; // *...*
	temp[11] = BH->charge*temp[4]; // *...*
	temp[12] = precalculated->rdot*chi_p; // *...*
	temp[13] = -p->mass*precalculated->phidot; // -...*
	temp[14] = p->charge*temp[11]; // +...*
	temp[15] = 0.5*BH->a5*p->charge; // +...*
	temp[16] = BH->a6*p->mass; // +...*
	temp[17] = p->r*precalculated->rdot; // *...*
	temp[18] = BH->mass*precalculated->phidot; // *...*
	temp[19] = temp[11]*sigma_p; // *...+
	temp[20] = p->mass*precalculated->phidot; // *...*
	temp[21] = -BH->a2*p->mass; // -...*
	temp[22] = p->r4*precalculated->rdot; // *...*
	temp[23] = temp[5]*precalculated->rdot; // *...*
	temp[24] = precalculated->tdot*sigma_p; // *...+
	temp[25] = temp[4]*precalculated->tdot; // *...*
	temp[26] = p->r3*precalculated->rdot; // *...*
	temp[27] = precalculated->tdot*sigma_p; // *...-
	temp[28] = -0.5*BH->a4*p->charge; // -...*
	temp[29] = -0.5*BH->a2*p->charge; // -...*
	temp[30] = -BH->a3*p->charge; // -...*
	temp[31] = precalculated->rdot*delta_p; // *...*
	temp[32] = temp[17]*delta_p; // *...*
	temp[33] = temp[18]*temp[4]; // *...*
	temp[34] = temp[11]*delta_p; // *...*
	temp[35] = BH->a2*temp[20]; // *...*
	temp[36] = -2.*temp[20]*temp[9]; // -...*
	temp[37] = 2.*BH->a3*p->mass; // +...*
	temp[38] = delta_p*sigma_p; // *...-
	temp[39] = BH->mass*temp[25]; // *...*
	temp[40] = p->mass*temp[26]; // *...*
	temp[41] = delta_p*sigma_p; // *...+
	temp[42] = 0.5*p->charge; // +...*
	temp[43] = chi_p*temp[41]; // *...+
	temp[44] = p->charge*BH->charge; // *...*
	temp[45] = temp[8]*precalculated->phidot; // +...*
	temp[46] = temp[20]*temp[26]; // *...*
	temp[47] = -BH->a*p->mass; // -...*
	temp[48] = precalculated->delta_p2*sigma_p; // *...-
	temp[49] = BH->a2+p->r2-delta_p; // (...)
	temp[50] = temp[0]*BH->charge; // +...*
	temp[51] = temp[1]*temp[5]; // +...*
	temp[52] = temp[2]*BH->charge; // +...*
	temp[53] = p->charge*temp[11]; // *...*
	temp[54] = p->charge*temp[23]; // *...*
	temp[55] = temp[11]*delta_p; // *...+
	temp[56] = temp[45]*temp[9]; // +...*
	temp[57] = precalculated->phidot*temp[10]; // *...*
	temp[58] = chi_p*delta_p; // *...-
	temp[59] = temp[5]*temp[12]; // *...*
	temp[60] = temp[11]*precalculated->delta_p2; // *...+
	temp[61] = temp[14]*chi_p; // +...*
	temp[62] = temp[15]*BH->charge; // +...*
	temp[63] = temp[16]*precalculated->phidot; // +...*
	temp[64] = p->mass*temp[33]; // *...*
	temp[65] = p->charge*temp[19]; // *...+
	temp[66] = temp[21]*temp[18]; // +...*
	temp[67] = 0.5*BH->a*temp[54]; // +...*
	temp[68] = temp[9]*sigma_p; // *...-
	temp[69] = temp[17]*temp[24]; // *...+
	temp[70] = temp[39]*sigma_p; // *...-
	temp[71] = temp[40]*temp[24]; // *...+
	temp[72] = BH->mass*temp[22]; // *...*
	temp[73] = temp[9]*temp[27]; // *...+
	temp[74] = temp[12]*sigma_p; // *...+
	temp[75] = chi_p*sigma_p; // *...+
	temp[76] = BH->charge*temp[31]; // *...*
	temp[77] = temp[20]*temp[32]; // *...*
	temp[78] = temp[33]*temp[41]; // *...+
	temp[79] = -5.*temp[35]*temp[26]; // -...*
	temp[80] = temp[37]*temp[17]; // +...*
	temp[81] = temp[47]*temp[39]; // +...*
	temp[82] = BH->a*temp[40]; // *...*
	temp[83] = BH->a2*temp[44]; // +...*
	temp[84] = temp[42]*temp[11]; // +...*
	temp[85] = BH->a*temp[44]; // *...*
	temp[86] = precalculated->delta_p2*sigma_p; // *...+
	temp[87] = 2.*temp[46]*temp[86]; // +...+
	temp[88] = precalculated->tdot*temp[48]; // *...-
	temp[89] = 0.+4.*p->r4*delta_p; // (...)
	temp[90] = temp[50]*p->r2; // +...*
	temp[91] = temp[51]*precalculated->rdot; // +...+
	temp[92] = temp[53]*chi_p; // *...+
	temp[93] = temp[6]*temp[55]; // +...+
	temp[94] = temp[56]*delta_p; // +...+
	temp[95] = -2.*BH->a2*temp[53]; // -...*
	temp[96] = temp[59]*delta_p; // *...+
	temp[97] = temp[13]*temp[9]; // +...*
	temp[98] = temp[61]*precalculated->delta_p2; // +...+
	temp[99] = temp[63]*temp[17]; // +...*
	temp[100] = temp[64]*sigma_p; // *...+
	temp[101] = 3.*BH->a4*temp[46]; // +...*
	temp[102] = temp[22]*sigma_p; // *...+
	temp[103] = 2.*temp[35]*temp[68]; // +...-
	temp[104] = BH->a3*p->mass; // +...*
	temp[105] = -2.*BH->a3*temp[71]; // -...+
	temp[106] = temp[72]*temp[27]; // *...+
	temp[107] = temp[28]*BH->charge; // +...*
	temp[108] = temp[11]*temp[75]; // *...+
	temp[109] = -2.*BH->a4*temp[77]; // -...*
	temp[110] = temp[2]*temp[34]; // +...*
	temp[111] = temp[36]*temp[41]; // +...+
	temp[112] = temp[81]*temp[41]; // +...+
	temp[113] = precalculated->tdot*temp[41]; // *...+
	temp[114] = temp[84]*temp[43]; // +...+
	temp[115] = precalculated->rdot*temp[86]; // *...+
	temp[116] = temp[87]+temp[47]*temp[17]*temp[88]; // +...-
	temp[117] = temp[49]*temp[89]; // *...*
	temp[118] = temp[90]*precalculated->rdot; // +...+
	temp[119] = BH->a4*temp[92]; // +...+
	temp[120] = temp[93]+temp[7]*temp[23]*delta_p; // +...+
	temp[121] = temp[95]*temp[58]; // +...-
	temp[122] = temp[2]*temp[60]; // +...+
	temp[123] = temp[98]+temp[62]*precalculated->rdot*sigma_p; // +...+
	temp[124] = temp[101]*sigma_p; // +...+
	temp[125] = temp[67]*sigma_p; // +...+
	temp[126] = -BH->a5*p->mass; // -...*
	temp[127] = temp[104]*temp[70]; // +...+
	temp[128] = p->mass*temp[106]; // *...+
	temp[129] = temp[107]*temp[74]; // +...+
	temp[130] = temp[30]*temp[76]; // +...*
	temp[131] = temp[8]*temp[78]; // +...+
	temp[132] = temp[79]*temp[41]; // +...+
	temp[133] = temp[80]*precalculated->tdot; // +...*
	temp[134] = temp[82]*temp[113]; // *...+
	temp[135] = temp[114]+0.5*temp[85]*temp[115]; // +...+
	temp[136] = -0.5*temp[44]*temp[12]; // -...*
	temp[137] = temp[118]+temp[91]; // +...+
	temp[138] = temp[119]+BH->a2*temp[54]*chi_p; // +...+
	temp[139] = temp[57]*delta_p; // *...+
	temp[140] = -p->charge*temp[96]; // -...+
	temp[141] = temp[123]+temp[99]*sigma_p; // +...-
	temp[142] = temp[124]+temp[66]*temp[102]; // +...+
	temp[143] = temp[127]+temp[105]; // +...+
	temp[144] = temp[47]*temp[73]; // +...+
	temp[145] = temp[29]*temp[108]; // +...+
	temp[146] = temp[109]*sigma_p; // +...+
	temp[147] = temp[110]*sigma_p; // +...+
	temp[148] = temp[111]+temp[133]*temp[38]; // +...+
	temp[149] = temp[83]*temp[12]; // +...*
	temp[150] = temp[45]*temp[17]; // +...*
	temp[151] = temp[136]*precalculated->delta_p2; // +...*
	temp[152] = temp[137]+temp[52]*temp[3]; // +...+
	temp[153] = temp[94]+p->mass*temp[139]; // +...+
	temp[154] = temp[122]+temp[97]*precalculated->delta_p2; // +...+
	temp[155] = temp[142]+temp[125]; // +...+
	temp[156] = temp[143]+BH->a*temp[128]; // +...+
	temp[157] = temp[145]+temp[130]*sigma_p; // +...+
	temp[158] = temp[147]+temp[132]; // +...+
	temp[159] = 2.*temp[134]+temp[149]*temp[41]; // +...+
	temp[160] = temp[152]+temp[138]; // +...+
	temp[161] = temp[121]+temp[140]; // +...+
	temp[162] = -BH->a4*temp[100]; // -...+
	temp[163] = temp[155]+temp[103]; // +...+
	temp[164] = temp[156]+temp[144]; // +...+
	temp[165] = temp[146]+temp[131]; // +...+
	temp[166] = temp[112]+temp[159]; // +...+
	temp[167] = temp[160]+temp[120]; // +...+
	temp[168] = temp[154]+temp[141]; // +...+
	temp[169] = temp[163]+temp[126]*temp[69]; // +...+
	temp[170] = temp[165]+temp[158]; // +...+
	temp[171] = temp[135]+temp[150]*temp[86]; // +...+
	temp[172] = temp[167]+temp[153]; // +...+
	temp[173] = temp[162]+BH->a3*temp[65]; // +...+
	temp[174] = temp[129]+temp[157]; // +...+
	temp[175] = temp[166]+temp[171]; // +...+
	temp[176] = temp[172]+temp[161]; // +...+
	temp[177] = temp[169]+temp[164]; // +...+
	temp[178] = temp[148]+temp[175]; // +...+
	temp[179] = temp[176]+temp[168]; // +...+
	temp[180] = temp[174]+temp[170]; // +...+
	temp[181] = temp[179]+temp[173]; // +...+
	temp[182] = temp[178]+temp[116]; // +...+
	temp[183] = temp[181]+temp[177]; // +...+
	temp[184] = temp[183]+temp[180]; // +...+
	temp[185] = temp[184]+temp[182]; // +...+
	temp[186] = +temp[185]+temp[151]*sigma_p; // (...)


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((8.*temp[186])/(p->mass*temp[117]*sigma_p));
	return phi_new;
}

cpp_dec_float_n calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n theta_new: New coordinate.
	*/

	cpp_dec_float_n theta_new{};


	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*(0.);
	return theta_new;
}

std::tuple<cpp_dec_float_n, cpp_dec_float_n,
	cpp_dec_float_n, cpp_dec_float_n> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n t_new: New t coordinate.
		- cpp_dec_float_n r_new: New r coordinate.
		- cpp_dec_float_n phi_new: New phi coordinate.
		- cpp_dec_float_n theta_new: New theta coordinate.
	*/

	cpp_dec_float_n t_new{};
	cpp_dec_float_n r_new{};
	cpp_dec_float_n phi_new{};
	cpp_dec_float_n theta_new{};

	Precalculated precalculated{BH, p, ps, sigma_p, delta_p, chi_p, dlambda};

	t_new = calc_t_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp0);

	r_new = calc_r_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp1);

	phi_new = calc_phi_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp2);

	theta_new = calc_theta_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	return {t_new, r_new, phi_new, theta_new};
}
