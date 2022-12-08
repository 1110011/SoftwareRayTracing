#ifndef LIGHTSCATTERINGDISTRIBUTIONSOLVER_H_INCLUDED
#define LIGHTSCATTERINGDISTRIBUTIONSOLVER_H_INCLUDED

// irradiance based on the Oren-Nayar reflectance model described in Oren et. al. "Generalization of Lambert's reflectance model" (1994)
double Li(vec3 &vIncidence, vec3 &rNormal, material mat, double theta_r, double azimuth)
{
    double cosTheta_i = dot(vIncidence,rNormal);
    double theta_i = acos(cosTheta_i);

    double alpha = std::max(theta_i, theta_r);
    double beta  = std::min(theta_i, theta_r);

    double A = 1-(0.5*(mat.sqRoughness/(mat.sqRoughness+0.33d)));
    double B = 0.45*(mat.sqRoughness/(mat.sqRoughness+0.09d));

    return (mat.albedo/M_PI) * cosTheta_i * (A + (B * std::max(0.0d, cos(azimuth)) * sin(alpha) * tan(beta)));
}

#endif // LIGHTSCATTERINGDISTRIBUTIONSOLVER_H_INCLUDED
