use std::fmt;

type Matrix = Vec<Vec<f64>>;
type Vector = Vec<f64>;

fn main() {
    let transition_matrix: Matrix = vec![
        vec![0.8, 0.1, 0.05], 
        vec![0.1, 0.7, 0.25], 
        vec![0.1, 0.2, 0.70], 
    ];

    let initial_vector: Vector = vec![1.0, 0.0, 0.0];

    println!("Simulation Started");
    println!("Initial State (t=0): {:.4?}", initial_vector);

    let mut current_vector = initial_vector;

    for step in 1..=5 {
        current_vector = multiply_matrix_vector(&transition_matrix, &current_vector);
        
        println!("State at Step {}:    {:.4?}", step, current_vector);
    }
    
    println!("Simulation Complete");
}

fn multiply_matrix_vector(matrix: &Matrix, vector: &Vector) -> Vector {
    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut result_vector = vec![0.0; rows];

    for i in 0..rows {
        let mut sum = 0.0;
        for j in 0..cols {
            sum += matrix[i][j] * vector[j];
        }
        result_vector[i] = sum;
    }

    result_vector
}

fn is_valid_dimension(matrix: &Matrix, vector: &Vector) -> bool {
    let rows = matrix.len();
    if rows == 0 { return false; }
    
    let cols = matrix[0].len();
    let vec_len = vector.len();

    rows == cols && cols == vec_len
}