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

import cc.redberry.groovy.feyncalc.pairedchi.SetupChi;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;

import java.util.List;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class RunChi implements Process {
    private final RunChiParameters parameters = new RunChiParameters();

    @Override
    public String command() {
        return "chi";
    }

    @Override
    public ProcessParameters parameters() {
        return parameters;
    }

    @Override
    public void run() {
        SetupChi.calc(parameters.getBottomSpin(), parameters.getCharmSpin(), parameters.getOutputFile());
    }

    @Parameters(commandDescription = "Run process g + g -> chi_b + chi_c",
            optionPrefixes = "-")
    public static final class RunChiParameters extends ProcessParameters {
        @Parameter(description = "bottomSpin charmSpin output_file")
        public List<String> params;

        @Override
        public String getOutputFileName() {
            return params.get(2);
        }

        public String getBottomSpin() {
            return params.get(0).toLowerCase();
        }

        public String getCharmSpin() {
            return params.get(1).toLowerCase();
        }

        @Override
        public void validate() {
            if (params.size() != 3)
                throw new ParameterException("Missing arguments.");
            for (String spin : new String[]{getBottomSpin(), getCharmSpin()})
                switch (spin) {
                    case "scalar":
                    case "axial":
                    case "tensor":
                        break;
                    default:
                        throw new ParameterException("Illegal spin argument: " + spin);
                }
            super.validate();
        }
    }
}
