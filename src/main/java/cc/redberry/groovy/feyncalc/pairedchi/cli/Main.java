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

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class Main {
    public static Process[] processes = {new RunCC(), new RunChi(), new RunCCXYZ(), new Convert()};
    public static Map<String, Process> processByName = new HashMap<>();

    static {
        for (Process p : processes)
            processByName.put(p.command(), p);
    }

    public static String mainCommand = "bbcc";

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            printHelp();
            return;
        }

        MainParameters mainParameters = new MainParameters();
        JCommander jCommander = new JCommander(mainParameters);
        jCommander.setProgramName(mainCommand);

        for (Process process : processes)
            jCommander.addCommand(process.command(), process.parameters());

        Process process = processByName.get(args[0]);
        if (process == null) {
            printHelp();
            return;
        }

        jCommander.parse(args);
        if (mainParameters.help()) {
            printHelp();
            return;
        }

        if (process.parameters().help()) {
            printProcessHelp(jCommander, process);
            return;
        }

        try {
            process.parameters().validate();
        } catch (ParameterException e) {
            System.out.println(e.getMessage());
            printProcessHelp(jCommander, process);
            return;
        }

        process.run();
    }

    public static void printHelp() {
        JCommander jCommander = new JCommander(new MainParameters());
        jCommander.setProgramName(mainCommand);
        for (Process process : processes)
            jCommander.addCommand(process.command(), process.parameters());
        StringBuilder builder = new StringBuilder();
        jCommander.usage(builder);
        System.out.println(builder);
    }


    public static void printProcessHelp(JCommander commander, Process process) {
        StringBuilder out = new StringBuilder();
        commander.usage(process.command(), out);
        System.out.println(out);
    }

    public static final class MainParameters {
        @Parameter(names = {"-h", "--help"}, help = true, description = "Displays this help message.")
        public Boolean help;

        public boolean help() {
            return help != null && help;
        }
    }
}
