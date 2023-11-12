/*
 * Copyright 2018-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package ProteomicsLibrary;

import ProteomicsLibrary.Types.Hypergeometric;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

import java.util.*;

public class Statistics {

    private static final NormalDistribution normalDistribution = new NormalDistribution(0, 1);

    public static double calMedian(Collection<Double> inputList) throws Exception {
        if (inputList.isEmpty()) {
            throw new Exception("There is no element in the input list.");
        } else {
            Double[] inputArray = inputList.toArray(new Double[0]);
            Arrays.sort(inputArray);
            if (inputArray.length % 2 == 0) {
                return (inputArray[inputArray.length / 2] + inputArray[(inputArray.length / 2) - 1]) / 2;
            } else {
                return inputArray[inputList.size() / 2];
            }
        }
    }

    public static double calMean(double[] input) throws Exception {
        if (input.length > 0) {
            double mean = 0;
            for (double v : input) {
                mean += v;
            }
            return mean / input.length;
        } else {
            throw new Exception("There is no element in the input array.");
        }
    }



    public static double calPearsonCorrelationCoefficient(double[] input1, double[] input2) throws Exception{
        if (input1.length != input2.length) {
            throw new Exception("Two vectors' lengths are different.");
        }

        double mean1 = calMean(input1);
        double mean2 = calMean(input2);
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        for (int i = 0; i < input1.length; ++i) {
            double c1 = input1[i] - mean1;
            double c2 = input2[i] - mean2;
            temp1 += c1 * c2;
            temp2 += Math.pow(c1, 2);
            temp3 += Math.pow(c2, 2);
        }
        return (temp1 == 0 || temp2 == 0) ? 0 : temp1 / (Math.sqrt(temp2 * temp3));
    }

}
