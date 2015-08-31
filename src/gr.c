/**
 * @file 	gr.c
 * @brief 	Post-newtonian general relativity corrections
 * @author 	Pengshuai (Sam) Shi <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "gr.h"

void rebx_gr(struct reb_simulation* const sim){
	struct rebx_params_gr* rebxparams = &((struct rebx_extras*)(sim->extras))->gr;
	const double C = rebxparams->c;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const particles = sim->particles;
	
	const struct reb_particle sun = particles[0];
	for (int i=1; i<_N_real; i++){
		const double dx = particles[i].x - sun.x;
		const double dy = particles[i].y - sun.y;
		const double dz = particles[i].z - sun.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double _r = sqrt(r2);
		const double dvx = particles[i].vx - sun.vx;
		const double dvy = particles[i].vy - sun.vy;
		const double dvz = particles[i].vz - sun.vz;
		// Benitez and Gallardo 2008
		const double alpha = G*sun.m/(_r*_r*_r*C*C);
		const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
		const double beta = 4.*G*sun.m/_r - v2;
		const double gamma = 4.*(dvx*dx + dvy*dy + dvz*dz);

		const double dax = alpha*(beta*dx + gamma*dvx);
		const double day = alpha*(beta*dy + gamma*dvy);
		const double daz = alpha*(beta*dz + gamma*dvz);
		const double massratio = particles[i].m/particles[0].m;

		particles[i].ax += dax;
		particles[i].ay += day;
		particles[i].az += daz;
		particles[0].ax -= massratio*dax;
		particles[0].ay -= massratio*day;
		particles[0].az -= massratio*daz;
	}
}

void rebx_gr_potential(struct reb_simulation* const sim){
	// Nobili & Roxburgh 1986
	struct rebx_params_gr* rebxparams = &((struct rebx_extras*)(sim->extras))->gr;
	const double C = rebxparams->c;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const particles = sim->particles;
	
	const struct reb_particle sun = particles[0];
	const double prefac1 = 6.*(G*sun.m)*(G*sun.m)/(C*C);
	for (int i=1; i<_N_real; i++){
		const double dx = particles[i].x - sun.x;
		const double dy = particles[i].y - sun.y;
		const double dz = particles[i].z - sun.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double prefac = prefac1/(r2*r2);
		
		particles[i].ax -= prefac*dx;
		particles[i].ay -= prefac*dy;
		particles[i].az -= prefac*dz;
	}
}

void rebx_gr_implicit(struct reb_simulation* const sim){
	struct rebx_params_gr* rebxparams = &((struct rebx_extras*)(sim->extras))->gr;
	const double C = rebxparams->c;
	const double C2i = 1./(C*C);
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const restrict particles = sim->particles;
	const unsigned int _gravity_ignore_10 = sim->gravity_ignore_10;

	if (rebxparams->allocatedN<_N_real){
		rebxparams->a_const = realloc(rebxparams->a_const, sizeof(struct reb_vec3d)*_N_real);
		rebxparams->a_newton = realloc(rebxparams->a_newton, sizeof(struct reb_vec3d)*_N_real);
		rebxparams->a_new = realloc(rebxparams->a_new, sizeof(struct reb_vec3d)*_N_real);
		rebxparams->a_old = realloc(rebxparams->a_old, sizeof(struct reb_vec3d)*_N_real);
		rebxparams->allocatedN = _N_real;
	}
	struct reb_vec3d* restrict const a_const = rebxparams->a_const;
	struct reb_vec3d* restrict const a_newton = rebxparams->a_newton;
	struct reb_vec3d* restrict a_new = rebxparams->a_new;
	struct reb_vec3d* restrict a_old = rebxparams->a_old;

	for (int i=0; i<_N_real; i++){
		a_newton[i].x = particles[i].ax;
		a_newton[i].y = particles[i].ay;
		a_newton[i].z = particles[i].az;
	}
	// extra newtonian terms
	if (_gravity_ignore_10 && _N_real>1){
		const double dx = particles[0].x - particles[1].x;
		const double dy = particles[0].y - particles[1].y;
		const double dz = particles[0].z - particles[1].z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double r = sqrt(r2);
		const double prefact = -G/(r2*r);
		const double prefact0 = prefact*particles[0].m;
		const double prefact1 = prefact*particles[1].m;
		a_newton[0].x = prefact1*dx;
		a_newton[0].y = prefact1*dy;
		a_newton[0].z = prefact1*dz;
		a_newton[1].x = prefact0*dx;
		a_newton[1].y = prefact0*dy;
		a_newton[1].z = prefact0*dz;

	}

	// then compute the constant terms:
	memset(a_const,0,sizeof(struct reb_vec3d)*_N_real);
	for (int i=0; i<_N_real; i++){
		// 1st constant part
		for (int j=0; j<_N_real; j++){
			if (j != i){
				double a1 = 0.;
				double a2 = 0.;
				for (int k=0; k<_N_real; k++){
					if (k != i){
						const double dxik = particles[i].x - particles[k].x;
						const double dyik = particles[i].y - particles[k].y;
						const double dzik = particles[i].z - particles[k].z;
						const double r2ik = dxik*dxik + dyik*dyik + dzik*dzik;
						const double rik = sqrt(r2ik);
						a1 += 4.* G*particles[k].m/rik;
					}
					if (k != j){
						const double dxlj = particles[k].x - particles[j].x;
						const double dylj = particles[k].y - particles[j].y;
						const double dzlj = particles[k].z - particles[j].z;
						const double r2lj = dxlj*dxlj + dylj*dylj + dzlj*dzlj;
						const double rlj = sqrt(r2lj);
						a2 += G*particles[k].m/rlj;
					}
				}
				
				const double dxij = particles[i].x - particles[j].x;
				const double dyij = particles[i].y - particles[j].y;
				const double dzij = particles[i].z - particles[j].z;
				const double r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
				const double rij = sqrt(r2ij);
				const double rij3i = 1./(r2ij*rij);
				
				const double dvxij = particles[i].vx - particles[j].vx;
				const double dvyij = particles[i].vy - particles[j].vy;
				const double dvzij = particles[i].vz - particles[j].vz;

				double vi2  = particles[i].vx*particles[i].vx 
				             +particles[i].vy*particles[i].vy 
					     +particles[i].vz*particles[i].vz;
				double a3 = -vi2;

				double vj2  = particles[j].vx*particles[j].vx 
				             +particles[j].vy*particles[j].vy 
					     +particles[j].vz*particles[j].vz;
				double a4 = -2.*vj2;

				double vij2 = particles[i].vx*particles[j].vx 
				             +particles[i].vy*particles[j].vy 
					     +particles[i].vz*particles[j].vz; 
				double a5 = 4.*vij2; 
				
				double dxijvj = dxij*particles[j].vx 
				               +dyij*particles[j].vy 
					       +dzij*particles[j].vz;
				double a6 = 3./2. * dxijvj*dxijvj/r2ij;
				
				double factor1 = a1 + a2 + a3 + a4 + a5 + a6;
					
				double factor2 = dxij*(4.*particles[i].vx-3.*particles[j].vx)
				                +dyij*(4.*particles[i].vy-3.*particles[j].vy)
						+dzij*(4.*particles[i].vz-3.*particles[j].vz);

				const double prefac1 = G*particles[j].m*factor1*rij3i*C2i;
				const double prefac2 = G*particles[j].m*factor2*rij3i*C2i;
				a_const[i].x += prefac1*dxij + prefac2*dvxij;
				a_const[i].y += prefac1*dyij + prefac2*dvyij;
				a_const[i].z += prefac1*dzij + prefac2*dvzij;
			}
		}
	}


	memset(a_new,0,sizeof(struct reb_vec3d)*_N_real);

	// Now running the iteration 
	for (int k=0; k<10; k++){ 
		{ // Swap
			struct reb_vec3d* restrict a_tmp = a_old;
			a_old = a_new;
			a_new = a_tmp;
			memset(a_new,0,sizeof(struct reb_vec3d)*_N_real);
		}
		// now add on the non-constant term
		for (int i=0; i<_N_real; i++){ // a_j is used to update a_i and vice versa
			for (int j=i+1; j<_N_real; j++){
				const double a_oldix = a_newton[i].x + a_const[i].x + a_old[i].x;
				const double a_oldiy = a_newton[i].y + a_const[i].y + a_old[i].y;
				const double a_oldiz = a_newton[i].z + a_const[i].z + a_old[i].z;
				                                        
				const double a_oldjx = a_newton[j].x + a_const[j].x + a_old[j].x;
				const double a_oldjy = a_newton[j].y + a_const[j].y + a_old[j].y;
				const double a_oldjz = a_newton[j].z + a_const[j].z + a_old[j].z;

				const double dxij = particles[i].x - particles[j].x;
				const double dyij = particles[i].y - particles[j].y;
				const double dzij = particles[i].z - particles[j].z;
				const double r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
				const double rij = sqrt(r2ij);
				const double daj = dxij*a_oldjx+dyij*a_oldjy+dzij*a_oldjz;
				const double dai = dxij*a_oldix+dyij*a_oldiy+dzij*a_oldiz;
				const double prefac1 = G/(r2ij*rij*2.*C*C);
				const double prefac2 = (7./2.*C2i)*G/rij;
				a_new[i].x += particles[j].m * (prefac1*daj*dxij + prefac2*a_oldjx);
				a_new[i].y += particles[j].m * (prefac1*daj*dyij + prefac2*a_oldjy);
				a_new[i].z += particles[j].m * (prefac1*daj*dzij + prefac2*a_oldjz);
				                                       
				a_new[j].x -= particles[i].m * (prefac1*dai*dxij + prefac2*a_oldix);
				a_new[j].y -= particles[i].m * (prefac1*dai*dyij + prefac2*a_oldiy);
				a_new[j].z -= particles[i].m * (prefac1*dai*dzij + prefac2*a_oldiz);
			}
		}

		// break out loop if a_new is converging
		double maxd = 0.;
		for (int i=0; i< _N_real; i++){
			const double dx = a_new[i].x - a_old[i].x;
			const double dy = a_new[i].y - a_old[i].y;
			const double dz = a_new[i].z - a_old[i].z;
			const double x = a_new[i].x;
			const double y = a_new[i].y;
			const double z = a_new[i].z;
			const double d = (dx*dx + dy*dy + dz*dz) / (x*x + y*y + z*z);
			if (isnormal(d)  && d>maxd){
				maxd = d;
			}
		}
		if (maxd<1e-30){
			break;
		}
	}
	// update acceleration in particles
	for (int i=0; i<_N_real;i++){
		particles[i].ax += a_new[i].x + a_const[i].x; 
		particles[i].ay += a_new[i].y + a_const[i].y;
		particles[i].az += a_new[i].z + a_const[i].z;
	}
					
}

