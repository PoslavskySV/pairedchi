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
package cc.redberry.groovy.feyncalc.pairedchi.cli;

import org.junit.Test;

/**
 * Created by poslavsky on 06/09/15.
 */
public class RunCCXYZTest {
    @Test
    public void test1() throws Exception {

        for (int i = 0; i < RunCCXYZ.encodedParams.length; i++) {
            RunCCXYZ.Encoded p = RunCCXYZ.encodedParams[i];
            System.out.println("xyz__" + i + ": " + p);
        }
    }
}