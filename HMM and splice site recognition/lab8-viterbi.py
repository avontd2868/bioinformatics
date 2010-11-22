def print_step(time_step,states,T):
   print '###############################################'
   print 'Time Step in Observation Sequence:', time_step
   for hidden_state in states:
      # T: (probability of hidden_state at current time, 
      #     [sequence of states that would lead to this hidden state, ie, Traceback], 
      #      highest probability contribution leading to this hidden state)
      print hidden_state, ': ', T[hidden_state] 

def viterbi(obs, states, start_p, trans_p, emit_p):
   verbose = False
   T = {}
   time_step = 0  # not needed for code, just for interpretation
   for state in states:
       ## initial transition probabilities
       ##          total prob , V-path(hidden state) , V-prob(hidden state)
       T[state] = (start_p[state], [state], start_p[state])
   if verbose == True:
      print_step(time_step,states,T)
   for output in obs:
       U = {}
       ## next_state is the hidden state at the current time point
       for next_state in states:  
           total = 0
           argmax = None
           valmax = 0
           ## source_state (hidden states at previous time point) 
           ## calculate the prob contributions from hidden states at previous time point 
           for source_state in states:  
               (prob, v_path, v_prob) = T[source_state]
               #print source_state, output, next_state
               #print emit_p[source_state][output]
               p = emit_p[source_state][output] * trans_p[source_state][next_state]
               prob *= p
               v_prob *= p
               total += prob
               #print("prob", prob)
               if v_prob > valmax:
                   argmax = v_path + [next_state]
                   valmax = v_prob
               #print("best prob: ",valmax)
               #print("best path: ",argmax)
           U[next_state] = (total, argmax, valmax)
       time_step = time_step + 1   # next observation in sequence
       T = U
       # debug, interpretation
       if verbose == True:
          print_step(time_step,states,T)       
   ## Traceback
   if verbose == True:
      print '###############################################'
      print '###############################################'
      print 'Printed are the possible Tracebacks from the final two possible hidden states.'
      print 'Code now picks the best path from the highest-probability final state (the first element of T list).'
      print_step(time_step,states,T)
   total = 0
   best_final_path = None
   valmax = 0
   print T
   for state in states:
       (state_prob, v_path, v_prob) = T[state]
       total += state_prob
       # start tracing back from the hidden state with highest final probability
       if state_prob > valmax:
           best_final_path = v_path
           valmax = state_prob
   return (total, best_final_path, valmax)
   # END Viterbi

def example_input1():
   states = ('E', '5p', 'I', '.')
   start_probability = {'E': 1.0, '5p':0, 'I':0, '.':0 }
   observations = ('C','T','T','C','A','T','G','T','G','A','A','A','G','C','A','G','A','C','G','T','A','A','G','T','C','A')
   transition_probability = {
      'E': {'E': 0.9, '5p': 0.1, 'I':0, '.':0 },
      '5p': {'E':0, '5p':0, 'I': 1.0, '.':0 },
      'I': {'E':0, '5p':0, 'I': 0.9, '.': 0.1 },
      '.': {'E':0, '5p':0, 'I': 0, '.': 0 },
   }
   emission_probability = {
      'E': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T':0.25 },
      '5p': {'A': 0.05, 'C': 0, 'G': 0.95, 'T':0 },
      'I': {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T':0.4 },
      '.': {'A': 0, 'C': 0, 'G': 0, 'T':0 },
      
   }
   return (states,start_probability,observations,transition_probability,emission_probability)

### Main ###
if __name__=='__main__':
   (states,start_probability,observations,transition_probability,emission_probability) = example_input1()
   (total_prob,final_path,final_state_prob) = viterbi(observations, states, start_probability, transition_probability, emission_probability)

   print 'Total probability of path ', total_prob, '.'
   print 'The final path, or sequence of hidden states:'
   print 'Observed:  ', observations
   print 'Predicted: ', final_path
   print 'Highest probabiilty final hidden state is ', final_path[-1], ' with relative probability ', final_state_prob, '.'
