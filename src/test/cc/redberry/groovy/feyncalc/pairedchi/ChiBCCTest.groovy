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

import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.groovy.RedberryPhysics.setMandelstam5

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCCTest {
    @Test
    public void testMc() throws Exception {
        def stp = new ChiBCC()
        println TensorUtils.getAllDiffSimpleTensors(stp.mc('scalar'))

    }

    @Test
    public void testCalculateAll() throws Exception {
        def stp = new ChiBCC()
        stp.calculateAll()
    }

    @Test
    public void test3() throws Exception {
        use(Redberry){
            def m = setMandelstam5([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]', 'p_m[bottom]': '2*m[bottom]'])

            println m
        }

    }
}
