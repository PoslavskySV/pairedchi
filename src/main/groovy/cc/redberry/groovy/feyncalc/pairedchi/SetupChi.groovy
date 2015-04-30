/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry

import static cc.redberry.groovy.RedberryStatic.Identity

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupChi extends Setup {
    private def qVertices
    private Map diagrams = [:]
    public Transformation epsVals = Identity

    SetupChi() {
        super(true, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values() as Transformation

            //polarizations of gluons
            def c1 = 'c1 = -(t-4*mb**2)*(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(-1/2)'.t
            def c2 = 'c2 = (t-4*mc**2+s)*(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(-1/2)'.t
            def c3 = 'c3 = -(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(1/2)'.t
            def c4 = 'c4 = -2*(-s*(t**2-4*t*mb**2+t*s+16*mc**2*mb**2-4*t*mc**2))**(-1/2)'.t
            def cc = c1 & c2 & c3 & c4
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t

            def epsPlus = 'epsPlus_a = (eps1_a + I * eps2_a)/2**(1/2)'.t
            def epsMinus = 'epsMinus_a = (eps1_a - I * eps2_a)/2**(1/2)'.t
            epsPlus <<= eps1 & eps2 & cc
            epsMinus <<= eps1 & eps2 & cc

            epsVals &= epsPlus >> 'eps1_a[+1] = epsPlus_a'.t
            epsVals &= epsMinus >> 'eps1_a[-1] = epsMinus_a'.t
            epsVals &= epsPlus >> 'eps2_a[+1] = epsPlus_a'.t
            epsVals &= epsMinus >> 'eps2_a[-1] = epsMinus_a'.t

            //polarizations of quarkonia (bottom)
            def s1 = 's1 = -4*mb*(16*mb**4-8*mb**2*s-32*mc**2*mb**2-8*mc**2*s+16*mc**4+s**2)**(-1/2)'.t
            def s2 = 's2 = -(1/2)*(4*mb**2+4*mc**2-s)*mb**(-1)*(16*mb**4-8*mb**2*s-32*mc**2*mb**2-8*mc**2*s+16*mc**4+s**2)**(-1/2)'.t
            def s3 = 's3 = (I)*(-4*mb**2*s-16*mc**2*mb**2-8*mc**2*s-4*t*mc**2+16*mc**4+s**2+4*t*mb**2+t*s)*((16*mc**2*mb**2-4*t*mc**2+t**2-4*t*mb**2+t*s)*s*(16*mb**4-8*mb**2*s-32*mc**2*mb**2-8*mc**2*s+16*mc**4+s**2))**(-1/2)'.t
            def s4 = 's4 = (-I)*(16*mc**2*mb**2+4*mc**2*s+4*t*mc**2-16*mc**4-4*t*mb**2+t*s)*((16*mc**2*mb**2-4*t*mc**2+t**2-4*t*mb**2+t*s)*s*(16*mb**4-8*mb**2*s-32*mc**2*mb**2-8*mc**2*s+16*mc**4+s**2))**(-1/2)'.t
            def s5 = 's5 = (-I)*(-4*mb**2-4*mc**2+8*mc*mb+s)*((16*mc**2*mb**2-4*t*mc**2+t**2-4*t*mb**2+t*s)*s*(16*mb**4-8*mb**2*s-32*mc**2*mb**2-8*mc**2*s+16*mc**4+s**2))**(-1/2)*(-4*mb**2-4*mc**2-8*mc*mb+s)'.t
            def s6 = 's6 = -2*(-(16*mc**2*mb**2-4*t*mc**2+t**2-4*t*mb**2+t*s)*s)**(-1/2)'.t
            def ss = s1 & s2 & s3 & s4 & s5 & s6

            eps1 = 'eps1_a = s1 * p_a[charm] + s2 * p_a[bottom]'.t
            eps2 = 'eps2_a = s3 * p_a[charm] + s4 * p_a[bottom] + s5 * k1_a'.t
            def eps0 = 'eps0_a = s6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]'.t
            epsPlus = 'epsPlus_a = (eps1_a + I * eps2_a)/2**(1/2)'.t
            epsMinus = 'epsMinus_a = (eps1_a - I * eps2_a)/2**(1/2)'.t
            epsPlus <<= eps1 & eps2; epsMinus <<= eps1 & eps2

            epsVals &= (epsPlus & ss) >> 'eps_a[bottom + 1] = epsPlus_a'.t
            epsVals &= (epsMinus & ss) >> 'eps_a[bottom - 1] = epsMinus_a'.t
            epsVals &= ((eps0 & ss) >> 'eps_a[bottom] = eps0_a'.t).hold

            //final tensor quarkonia

            //polarizations of quarkonia (charm)
            s1 = 's1 = -(1/2)*mc**(-1)*(-s+4*mc**2+4*mb**2)*(-8*s*mc**2-8*s*mb**2+s**2+16*mc**4-32*mc**2*mb**2+16*mb**4)**(-1/2)'.t
            s2 = 's2 = -4*mc*(-8*s*mc**2-8*s*mb**2+s**2+16*mc**4-32*mc**2*mb**2+16*mb**4)**(-1/2)'.t
            s3 = 's3 = (I)*(-8*s*mc**2+s*t-4*s*mb**2+s**2+16*mc**4-4*mc**2*t-16*mc**2*mb**2+4*t*mb**2)*(s*(t**2+s*t-4*mc**2*t+16*mc**2*mb**2-4*t*mb**2)*(-8*s*mc**2-8*s*mb**2+s**2+16*mc**4-32*mc**2*mb**2+16*mb**4))**(-1/2)'.t
            s4 = 's4 = (-I)*(4*s*mc**2+s*t-16*mc**4+4*mc**2*t+16*mc**2*mb**2-4*t*mb**2)*(s*(t**2+s*t-4*mc**2*t+16*mc**2*mb**2-4*t*mb**2)*(-8*s*mc**2-8*s*mb**2+s**2+16*mc**4-32*mc**2*mb**2+16*mb**4))**(-1/2)'.t
            s5 = 's5 = (-I)*(s-4*mc**2-4*mb**2+8*mb*mc)*(s-4*mc**2-4*mb**2-8*mb*mc)*(s*(t**2+s*t-4*mc**2*t+16*mc**2*mb**2-4*t*mb**2)*(-8*s*mc**2-8*s*mb**2+s**2+16*mc**4-32*mc**2*mb**2+16*mb**4))**(-1/2)'.t
            s6 = 's6 = -2*(-s*(t**2+s*t-4*mc**2*t+16*mc**2*mb**2-4*t*mb**2))**(-1/2)'.t
            ss = s1 & s2 & s3 & s4 & s5 & s6
            epsVals &= (epsPlus & ss) >> 'eps_a[charm + 1] = epsPlus_a'.t
            epsVals &= (epsMinus & ss) >> 'eps_a[charm - 1] = epsMinus_a'.t
            epsVals &= ((eps0 & ss) >> 'eps_a[charm] = eps0_a'.t).hold

            this.fullSimplify &= 'e_abcd * k1^a * k2^b * p^c[charm] * p^d[bottom] = 0'.t
        }
    }

    Tensor setupFeynmanDiagrams(charmSpin, bottomSpin) {
        if (diagrams[charmSpin + bottomSpin] != null)
            return diagrams[charmSpin + bottomSpin]

        use(Redberry) {
            log "Setting up diagrams for charm: $charmSpin bottom: $bottomSpin ..."
            def simplify = FeynmanRules & qVertices & fullSimplify
            //first diagram
            def Ma = '1'.t
            Ma *= simplify >> "eps1^a[h1] * A${charmSpin}_{aA cC}[charm, k1_i, -k1_i + p_i[charm]]".t
            Ma *= simplify >> 'G^cd[k1_i - p_i[charm]] * g^CD'.t
            Ma *= simplify >> "eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t

            //second diagram
            def Mb = '1'.t
            Mb *= simplify >> "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]]".t
            Mb *= simplify >> 'G^cd[k1_i - p_i[bottom]] * g^CD'.t
            Mb *= simplify >> "eps2^b[h2] * A${charmSpin}_{dD bB}[charm, -k2_i + p_i[charm], k2_i]".t

            //total matrix element
            def M = Ma + Mb
            M <<= momentumConservation
            M <<= fullSimplify & massesSubs & mFactor

            log '...done'
            return (diagrams[charmSpin + bottomSpin] = M)
        }
    }
//
//    public static void calc(String bottomSpin, String charmSpin, File output) {
//        use(Redberry) {
//            SetupChi stp = new SetupChi();
//            def amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
//            def amp2 = stp.squareMatrixElement(amps)
//            amp2 <<= stp.mandelstam & stp.massesSubs
//            def r = "charm${charmSpin}bottom${bottomSpin}".t.eq(amp2)
//            output << r.toString(OutputFormat.Redberry)
//        }
//    }

    public static void calc(def g1, def g2, String bottomSpin, String charmSpin, File output) {
        use(Redberry) {
            checkPol g1
            checkPol g2

            SetupChi stp = new SetupChi();
            def amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
            if (g1 != null)
                amps <<= "h1 = $g1".t
            if (g1 != null)
                amps <<= "h2 = $g2".t
            if (g1 != null || g2 != null) {
                stp.log "Setting up gluon polarizations $g1 $g2 ..."
                amps <<= stp.epsVals & stp.fullSimplifyE & stp.mFactor
                stp.log '... done'
            }

            def amp2 = stp.squareMatrixElement(amps)

            amp2 <<= stp.mandelstam & stp.massesSubs
            def r = "charm${charmSpin}bottom${bottomSpin}".t.eq(amp2)
            output << r.toString(OutputFormat.Redberry)
        }
    }

    private static void checkPol(g) {
        if (g != null && g != 1 && g != -1)
            throw new IllegalArgumentException()
    }
}
