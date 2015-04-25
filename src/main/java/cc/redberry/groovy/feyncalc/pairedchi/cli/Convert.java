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
import cc.redberry.core.tensor.Tensors;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;

import java.io.BufferedReader;
import java.io.InputStreamReader;
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
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(System.in));
        String line;
        OutputFormat of = parameters.getOutputFormat();
        while ((line = bufferedReader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty())
                continue;
            if (line.endsWith(";"))
                line = line.substring(0, line.length() - 1);
            try {
                System.out.println(Tensors.parse(line).toString(of));
            } catch (Exception e) {
                System.err.println("Invalid input: " + line);
                System.err.println(e.getMessage());
                System.exit(1);
            }
        }
    }

    @Parameters(commandDescription = "Converts Redberry output to Mathematica or Maple",
            optionPrefixes = "-")
    public static final class ConvertParameters extends ProcessParameters {
        @Parameter(description = "format input_stream")
        public List<String> params;

        @Override
        public String getOutputFileName() {
            return null;
        }

        public OutputFormat getOutputFormat() {
            switch (params.get(0).toLowerCase()) {
                case "w":
                case "math":
                case "mathematica":
                    return OutputFormat.WolframMathematica;
                case "m":
                case "maple":
                    return OutputFormat.Maple;
                default:
                    throw new ParameterException("Wrong format");
            }
        }

        @Override
        public void validate() {
        }
    }
}
