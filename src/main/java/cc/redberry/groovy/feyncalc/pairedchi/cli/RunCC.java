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

import java.util.List;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class RunCC implements Process {
    private final RunCCParameters parameters = new RunCCParameters();

    @Override
    public String command() {
        return "cc";
    }

    @Override
    public ProcessParameters parameters() {
        return parameters;
    }

    @Override
    public void run() {
        if (parameters.ward())
            SetupCC.ward(parameters.getSpin(), parameters.getOutputFile());
        else
            SetupCC.calc(parameters.getSpin(), parameters.getOutputFile());
    }

    @Parameters(commandDescription = "Run process g + g -> chi_b + c + cBar",
            optionPrefixes = "-")
    public static final class RunCCParameters extends ProcessParameters {
        @Parameter(description = "bottomSpin output_file")
        public List<String> params;

        @Parameter(names = {"-w", "--ward"}, description = "Ward check only: replace gluon polarizations with momentums")
        public Boolean ward;

        @Override
        public String getOutputFileName() {
            return params.get(1);
        }

        public String getSpin() {
            return params.get(0).toLowerCase();
        }

        public boolean ward() {
            return ward != null && ward;
        }

        @Override
        public void validate() {
            if (params.size() != 2)
                throw new ParameterException("Missing arguments.");
            switch (getSpin()) {
                case "scalar":
                case "axial":
                case "tensor":
                    break;
                default:
                    throw new ParameterException("Illegal spin argument: " + getSpin());
            }
            super.validate();
        }
    }
}
