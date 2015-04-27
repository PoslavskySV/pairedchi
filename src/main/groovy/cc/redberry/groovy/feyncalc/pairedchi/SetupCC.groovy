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

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry

import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class SetupCC extends Setup {
    private def qVertices, ccVertex

    private Map quarkDiagrams = [:]
    private Map gluonDiagrams = [:]
    private Map gluon3Diagrams = [:]
    private Tensor quarkDiagramsNotProjected = null;
    private Tensor gcc = null;

    private Transformation spinorsSimplify, dSimplify
    private Transformation assumptions

    SetupCC() {
        this(Identity)
    }

    SetupCC(assumptions) {
        super(false, true);
        this.assumptions = assumptions
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values() as Transformation
            ccVertex = effectivePairVertex()
            spinorsSimplify = Identity
            spinorsSimplify &= 'cu[p1_m[charm]]*G_a*p1^a[charm] = m[charm]*cu[p1_m[charm]]'.t
            spinorsSimplify &= 'G_a*p2^a[charm]*v[p2_m[charm]] = -m[charm]*v[p2_m[charm]]'.t
            dSimplify = Identity
            dSimplify &= 'G_a*G^a = 4'.t
            dSimplify &= 'G_a*G_b*G^a = -2*G_b'.t
            dSimplify &= 'G_a*G_b*G_c*G^a = 4*g_bc'.t

            gcc = 'Vcc_iI = G_ij[p1_m[charm] + p2_m[charm]]*cu[p1_m[charm]]*V^j_I*v[p2_m[charm]]'.t
            gcc <<= FeynmanRules & fullSimplify
        }
    }

    Tensor getGluonDiagrams(bottomSpin) {
        if (gluonDiagrams[bottomSpin] != null)
            return gluonDiagrams[bottomSpin]
        use(Redberry) {
            log "Setting up gluon diagrams for $bottomSpin ..."
            def simplify = assumptions & FeynmanRules & qVertices & ccVertex & fullSimplify
            //first diagram
            def Ma = 1.t
            Ma *= simplify >> 'eps1^a[h1] * B_{aA cC}[charm, k1_i, k2_i - p_i[bottom]]'.t
            Ma *= simplify >> 'G^cd[-k2_i + p_i[bottom]] * g^CD'.t
            Ma *= simplify >> "eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t

            //second diagram
            def Mb = 1.t
            Mb *= simplify >> 'eps2^b[h2] * B_{bB cC}[charm, k2_i, k1_i - p_i[bottom]]'.t
            Mb *= simplify >> 'G^cd[-k1_i + p_i[bottom]] * g^CD'.t
            Mb *= simplify >> "eps1^a[h1] * A${bottomSpin}_{dD aA}[bottom, -k1_i + p_i[bottom], k1_i]".t

            def M = Ma + Mb
            M <<= fullSimplify & massesSubs & mFactor & spinorsSimplify & massesSubs & mFactor

            log "... done"
            return (gluonDiagrams[bottomSpin] = M)
        }
    }

    Tensor get3GluonDiagrams(bottomSpin) {
        if (gluon3Diagrams[bottomSpin] != null)
            return gluon3Diagrams[bottomSpin]
        use(Redberry) {
            log "Setting up 3-gluon diagrams for $bottomSpin ..."
            def simplify = assumptions & FeynmanRules & qVertices & ccVertex & fullSimplify

            //first diagram
            log 'a) ...'
            def Ma = '1'.t
            Ma *= simplify >> 'eps1^a[h1] * Vcc^cC * V_{cC aA dD}[-p1_a[charm] - p2_a[charm], k1_a, k2_a - p_a[bottom]]'.t
            Ma *= simplify >> "G^de[-k2_a + p_a[bottom]]*g^DE".t
            Ma *= simplify >> "A${bottomSpin}_{eE bB}[bottom, -k2_a + p_a[bottom], k2_a] * eps2^b[h2]".t

            //second diagram
            log 'b) ...'
            def Mb = '1'.t
            Mb *= simplify >> "eps2^b[h2] * Vcc^cC * V_{cC bB dD}[-p1_a[charm] - p2_a[charm], k2_a, k1_a - p_a[bottom]]".t
            Mb *= simplify >> "G^ed[-k1_a + p_a[bottom]]*g^ED".t
            Mb *= simplify >> "A${bottomSpin}_{eE aA}[bottom, -k1_a + p_a[bottom], k1_a] * eps1^a[h1]".t

            //third diagram (s-chanel)
            log 'c) ...'
            def Mc = '1'.t
            Mc *= simplify >> 'eps1^a[h1] * eps2^b[h2] * V_{aA bB cC}[k1_a, k2_a, -k1_a - k2_a]'.t
            Mc *= simplify >> "G^cd[k1_a + k2_a]*g^CD".t
            Mc *= simplify >> "A${bottomSpin}_{dD eE}[bottom, k1_a + k2_a, -p1_a[charm] - p2_a[charm]] * Vcc^eE".t


            log 'simplifying ...'
            def M = Ma + Mb + Mc
            M <<= fullSimplify & massesSubs & mFactor & spinorsSimplify & massesSubs & mFactor
            M <<= (gcc << (massesSubs & mFactor))
            log 'simplifying ...'
            M <<= fullSimplifyE & massesSubs & mFactor & spinorsSimplify & massesSubs & mFactor

            log "... done"
            return (gluon3Diagrams[bottomSpin] = M)
        }
    }

    Tensor getQuarkDiagramsNotProjected() {
        if (quarkDiagramsNotProjected != null)
            return quarkDiagramsNotProjected

        use(Redberry) {

            quarkDiagramsNotProjected = 0.t
            // (1,2,3)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_bB*eps2^b[h2]*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,2,1)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_bB*eps2^b[h2]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (1,3,2)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_cC*Vcc^cC*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,1,2)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_aA*eps1^a[h1]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (2,3,1)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_cC*Vcc^cC*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t
            // (2,1,3)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_aA*eps1^a[h1]*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t


            def masses = 'p2_{f}[bottom]*p2^{f}[bottom] = m[bottom]**2'.t & 'p1_{d}[bottom]*p1^{d}[bottom] = m[bottom]**2'.t
            quarkDiagramsNotProjected <<= assumptions
            quarkDiagramsNotProjected <<= 'pCharm_m = p1_m[charm] + p2_m[charm]'.t
            quarkDiagramsNotProjected <<= FeynmanRules & spinSingletProjector['bottom']
            quarkDiagramsNotProjected <<= dTraceSimplify
            quarkDiagramsNotProjected <<= masses
            quarkDiagramsNotProjected <<= momentums['bottom']
            quarkDiagramsNotProjected <<= 'q_i[bottom] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[bottom]'.t
            quarkDiagramsNotProjected <<= mandelstam & ExpandAll[EliminateMetrics & mandelstam] & uTrace & EliminateMetrics
            return quarkDiagramsNotProjected
        }
    }

    Tensor getQuarkDiagrams(String bottomSpin) {
        if (quarkDiagrams[bottomSpin] != null)
            return quarkDiagrams[bottomSpin]
        use(Redberry) {
            log "Setting up quark diagrams for $bottomSpin ..."
            def Mc = getQuarkDiagramsNotProjected()
            Mc <<= totalSpinProjector[bottomSpin]
            if (bottomSpin == 'axial')
                Mc <<= momentumConservation
            log 'simplifying ...'
            Mc <<= fullSimplify & massesSubs & mFactor
            Mc <<= (gcc << (massesSubs & mFactor))
            log 'simplifying ...'
            Mc <<= fullSimplifyE & massesSubs & mFactor
            Mc <<= spinorsSimplify & massesSubs & mFactor
            log "... done"
            return (quarkDiagrams[bottomSpin] = Mc)
        }
    }

    public static void ward(String bottomSpin, File output) {
        for (int k = 1; k <= 2; ++k) {
            CC.reset()
            println "k = $k"

            use(Redberry) {
                def epss = "eps${k}_a[h${k}] = k${k}_a".t
                SetupCC stp = new SetupCC(epss);
                def methods = [1: stp.&getGluonDiagrams, 0: stp.&getQuarkDiagrams, 2: stp.&get3GluonDiagrams]
                def amps = new SumBuilder()

                for (int j = 0; j < 3; ++j) {
                    def diag = methods[j](bottomSpin)
                    diag <<= (epss & stp.mandelstam & stp.massesSubs)
                    amps << diag
                }
                def amp2 = stp.squareMatrixElement(amps.build())
                amp2 <<= stp.mandelstam & stp.massesSubs
                def r = "ward$bottomSpin$k".t.eq(amp2)
                output << r.toString(OutputFormat.Redberry)
                output << '\n'
            }
        }
    }

    public static void calc(String bottomSpin, File output) {
        use(Redberry) {
            SetupCC stp = new SetupCC();
            def methods = [1: stp.&getGluonDiagrams, 0: stp.&getQuarkDiagrams, 2: stp.&get3GluonDiagrams]
            def amps = new SumBuilder()

            for (def i = 0; i < 3; ++i)
                amps << methods[i](bottomSpin)

            def amp2 = stp.squareMatrixElement(amps.build())
            amp2 <<= stp.mandelstam & stp.massesSubs
            def r = "r$i".t.eq(amp2)
            output << r.toString(OutputFormat.Redberry)
        }
    }
}
