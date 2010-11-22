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
   states = ('Rainy', 'Sunny')
   start_probability = {'Rainy': 0.6, 'Sunny': 0.4}
   observations = ('walk', 'shop', 'clean', 'walk')
   transition_probability = {
   'Rainy' : {'Rainy': 0.7, 'Sunny': 0.3},
   'Sunny' : {'Rainy': 0.4, 'Sunny': 0.6},
   }
   emission_probability = {
   'Rainy' : {'walk': 0.1, 'shop': 0.4, 'clean': 0.5},
   'Sunny' : {'walk': 0.6, 'shop': 0.3, 'clean': 0.1},
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
