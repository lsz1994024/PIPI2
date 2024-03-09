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

package ProteomicsLibrary.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Double> sparseVector = new HashMap<>();

    public SparseVector() {}
    public void add(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            if (sparseVector.containsKey(i)) {
                sparseVector.put(i, sparseVector.get(i) + v);
            } else {
                sparseVector.put(i, v);
            }
        }
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public double get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }
    public Double[] getValues() {
        return sparseVector.values().toArray(new Double[0]);
    }

    public boolean isEmpty() {
        return sparseVector.isEmpty();
    }
}
