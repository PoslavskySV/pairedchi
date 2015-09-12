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

import cc.redberry.groovy.feyncalc.pairedchi.SetupCC;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class RunCCXYZ implements Process {
    private final RunCCXYZParameters parameters = new RunCCXYZParameters();

    @Override
    public String command() {
        return "xyz";
    }

    @Override
    public ProcessParameters parameters() {
        return parameters;
    }

    @Override
    public void run() {
        SetupCC stp = new SetupCC(true);
        int[] codes = parameters.getCodes();
        Map map = new HashMap();
        for (int code : codes) {
            Encoded p = encodedParams[code];
            File output = new File(parameters.getOutputFileName() + "__" + code);
            output.delete();
            SetupCC.calcXYZ(stp, p.g1, p.g2, p.S, p.L, output, map);
        }
    }

    static final class Encoded {
        final int g1, g2;
        final String S, L;

        public Encoded(int g1, int g2, String s, String l) {
            this.g1 = g1;
            this.g2 = g2;
            S = s;
            L = l;
        }


        @Override
        public String toString() {
            return "g1=" + g1 +
                    " g2=" + g2 +
                    " S=" + S +
                    " L=" + L;
        }
    }

    static final Encoded[] encodedParams;

    static {
        encodedParams = new Encoded[36];
        int[] gs = {-1, 1};
        String[] xyzs = {"x", "y", "z"};
        int i = 0;
        for (int g1 : gs)
            for (int g2 : gs)
                for (String s : xyzs)
                    for (String l : xyzs)
                        encodedParams[i++] = new Encoded(g1, g2, s, l);
    }

    @Parameters(commandDescription = "Calc amplitudes for g + g -> chi_b + c + cBar in XYZ plane",
            optionPrefixes = "-")
    public static final class RunCCXYZParameters extends ProcessParameters {
        @Parameter(description = "output_file")
        public List<String> params;

        @Parameter(names = {"-c", "--codes"}, description = "Parameters set")
        public List<String> codes;

        @Override
        public String getOutputFileName() {
            return params.get(0);
        }

        public int[] getCodes() {
            int[] codes = new int[this.codes.size()];
            for (int i = 0; i < codes.length; ++i)
                codes[i] = Integer.valueOf(this.codes.get(i));
            return codes;
        }

        @Override
        public void validate() {
            if (params.size() != 1)
                throw new ParameterException("No output file specified.");
            for (String codeString : codes) {
                int code = Integer.valueOf(codeString);
                if (code < 0 || code >= 36)
                    throw new ParameterException("Each code must be from 0 to 36.");
            }
            super.validate();
        }
    }
}

