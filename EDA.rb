class PBIL

    def initialize(n, iter_count, domains)
      @n = n
      @iter_count = iter_count
      @domains = domains
      #@costFunc = costFunc
      @learnRate = 0.1
      @negLearnRate = 0.075
      @mutProb = 0.02
      @mutShift = 0.05
      @totalBits = domains.size

      @bestCost
      @bestCosts = Array.new(@iter_count){Float::INFINITY}
      @bestGene = nil
    end


    def optimize
      @bestCost = Float::INFINITY
      probVec = Array.new(@totalBits){0.5}

      @iter_count.times do |i|

        # 個体生成
        genes = Array.new(@n){Array.new(@totalBits)}

        genes.each do |gene|
          gene.length.times do |k|
            gene[k] = (rand < probVec[k]) ? 1 : 0
          end
        end

        # コスト計算
        costs = Array.new(@n)
        genes.each_with_index do |gene,idx|
          costs[idx] = costFunc(toRealValue(gene))
        end

        # コスト最小・最大の個体を探す
        minGene = nil
        maxGene = nil
        minCost = Float::INFINITY
        maxCost = -Float::INFINITY
        @n.times do |j|
          cost = costs[j]
          if minCost > cost
            minCost = cost
            minGene = genes[j]
          end
          if maxCost < cost
            maxCost = cost
            maxGene = genes[j]
          end
        end

        # 全体の最小と比較
        if @bestCost > minCost
          @bestCost = minCost
          @bestGene = minGene
        end
        @bestCosts[i] = @bestCost

#        puts "mincost : #{minCost} minGene: #{minGene}"

        # 最小・最大コストで確率ベクトルを更新
        @totalBits.times do |j|
          if minGene[j] == maxGene[j]
            probVec[j] = probVec[j] * (1.0 - @learnRate) + (minGene[j] * @learnRate)
          else
            learnRate2 = @learnRate + @negLearnRate
            probVec[j] = probVec[j] * (1.0 - learnRate2) + (minGene[j] * learnRate2)
          end
        end

        # 突然変異
        @totalBits.times do |j|
          if rand < @mutProb
            probVec[j] = probVec[j] * (1.0 - @mutShift) + (rand(1) * @mutShift)
          end
        end

        p @bestCosts
      end
    end

  def toRealValue(gene)
    realValue = 0.0
    gene.each_with_index do |bit, idx|
      if idx == 0
        realValue += bit
      else
        realValue += 2**(idx*-1) * bit
      end
    end

    return realValue
  end

  def costFunc(x)
    0.993851231 + (Math::E ** (-0.01 * (x ** 2)) ) * Math.sin(10 * x) * Math.cos(8 * x)
  end

end

# -1.0～1.0まで
# 解像度2^24 -> -1.000000000000000000000000
# Array.new(25) ? -> 整数部[0/1] 実数部[0/1]*24
# -> Array.new(2){Array.new(25)} ... x と y
# probVec[25]
# 0.123

pbil = PBIL.new(100,100,Array.new(25))

pbil.optimize

p "hell"