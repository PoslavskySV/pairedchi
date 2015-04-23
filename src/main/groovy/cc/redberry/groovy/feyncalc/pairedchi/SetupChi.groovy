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

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupChi extends Setup {
    private def qVertices
    private Map diagrams = [:]

    SetupChi() {
        super(true, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values() as Transformation
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


    public static void calc(String bottomSpin, String charmSpin, File output) {
        use(Redberry) {
            SetupChi stp = new SetupChi();
            def amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
            def amp2 = stp.squareMatrixElement(amps)
            amp2 <<= stp.mandelstam & stp.massesSubs
            def r = "charm${charmSpin}bottom${bottomSpin}".t.eq(amp2)
            output << r.toString(OutputFormat.Redberry)
        }
    }
}
