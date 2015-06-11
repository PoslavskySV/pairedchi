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
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry

import static cc.redberry.groovy.RedberryStatic.Together

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupChi extends Setup {
    def qVertices

    SetupChi() {
        super(true, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values().collect({ x -> x << (Together & wFactor) }) as Transformation
        }
    }

    def diagrams(charmSpin, bottomSpin) {
        use(Redberry) {
            def Ma = "eps1^a[h1] * A${charmSpin}_{aA cC}[charm, k1_i, -k1_i + p_i[charm]] * G^cd[k1_i - p_i[charm]] * g^CD * eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t
            def Mb = "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]] * G^cd[k1_i - p_i[bottom]] * g^CD * eps2^b[h2] * A${charmSpin}_{dD bB}[charm, -k2_i + p_i[charm], k2_i]".t
            return [Ma, Mb]
        }
    }

    public static void calc(charmSpin, bottomSpin, g1, g2, File output) {
        use(Redberry) {
            SetupChi stp = new SetupChi();
            def diags = stp.diagrams(charmSpin, bottomSpin)
            def pols = stp.setupPolarisations(g1, g2)
            def M2 = stp.calcProcess(diags, pols)
            output << M2.toString(OutputFormat.Redberry)
        }
    }
}
