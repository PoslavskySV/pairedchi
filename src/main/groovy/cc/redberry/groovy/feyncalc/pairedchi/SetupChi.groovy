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
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry

import static cc.redberry.groovy.RedberryStatic.Identity
import static cc.redberry.groovy.RedberryStatic.Together

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupChi extends Setup {
    def qVertices
    private Map diagrams = [:]

    SetupChi() {
        super(true, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values().collect({ x -> x << (Together & mFactor) }) as Transformation
        }
    }

    Tensor setupFeynmanDiagrams(charmSpin, bottomSpin, apply = Identity) {
        if (diagrams[charmSpin + bottomSpin] != null)
            return diagrams[charmSpin + bottomSpin]

        use(Redberry) {
            log "Setting up diagrams for charm: $charmSpin bottom: $bottomSpin ..."
            def simplify = FeynmanRules & qVertices & apply & fullSimplify
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

            Ma <<= simplify
            Mb <<= simplify

            //total matrix element
            def M = Ma + Mb
            M <<= momentumConservation
            M <<= massesSubs & mFactor

            log '...done'
            return (diagrams[charmSpin + bottomSpin] = M)
        }
    }

    public static void calc(def g1, def g2, def charmPol, def bottomPol,
                            String bottomSpin, String charmSpin,
                            File output) {
        use(Redberry) {

            SetupChi stp = new SetupChi();
            def pols = stp.setPolarizations(g1, g2, charmPol, bottomPol, bottomSpin, charmSpin)
            def amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin, pols)
            println amps.class
            println amps.size()
            def sss = new SumBuilder()
            def c = 0
            for (def t in amps) {

                println t
                t *= 'g^AB'.t
                t <<= pols
                t <<= stp.fullSimplifyE & stp.massesSubs & stp.mFactor
                sss << t

                println(c++)
            }

            if (TensorUtils.isSymbolic(amps)) {
                stp.log 'scalar amplitude'
                amps <<= stp.massesSubs
                def r = "charm${charmSpin}bottom${bottomSpin}".t.eq(amps)
                output << r.toString(OutputFormat.Redberry)
                return
            } else {
                if (amps.class == Product) {
                    println TensorUtils.isSymbolic(t.indexlessSubProduct)
                    println t.dataSubProduct
                } else
                    for (def t in amps) {
                        if (t.class == Product)
                            println t.dataSubProduct
                    }
            }

            def amp2 = stp.squareMatrixElement(amps)
            amp2 <<= stp.mandelstam & stp.massesSubs
            def r = "charm${charmSpin}bottom${bottomSpin}".t.eq(amp2)
            output << r.toString(OutputFormat.Redberry)
        }
    }
}
