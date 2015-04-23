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

import cc.redberry.core.context.OutputFormat;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.util.List;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class Convert implements Process {
    private final ConvertParameters parameters = new ConvertParameters();

    @Override
    public String command() {
        return "convert";
    }

    @Override
    public ProcessParameters parameters() {
        return parameters;
    }

    @Override
    public void run() throws Exception {
        String redberry = FileUtils.readFileToString(new File(parameters.getInputFile()));
        Tensor t = Tensors.parse(redberry);
        try (FileOutputStream out = new FileOutputStream(parameters.getOutputFile())) {
            out.write(t.toString(parameters.getOutputFormat()).getBytes());
        }
    }

    @Parameters(commandDescription = "Converts Redberry output to Mathematica or Maple",
            optionPrefixes = "-")
    public static final class ConvertParameters extends ProcessParameters {
        @Parameter(description = "format input output")
        public List<String> params;

        public String getInputFile() {
            return params.get(1);
        }

        @Override
        public String getOutputFileName() {
            return params.get(2);
        }

        public OutputFormat getOutputFormat() {
            switch (params.get(0).toLowerCase()) {
                case "mathematica":
                    return OutputFormat.WolframMathematica;
                case "maple":
                    return OutputFormat.Maple;
                default:
                    throw new ParameterException("Wrong format");
            }
        }

        public String getCharmSpin() {
            return params.get(1).toLowerCase();
        }

        @Override
        public void validate() {
            if (params.size() != 3)
                throw new ParameterException("Missing arguments.");
            if (!new File(getInputFile()).exists())
                throw new ParameterException("Input file not exists.");
            getOutputFormat();
            super.validate();
        }
    }
}
