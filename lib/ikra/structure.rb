require "ikra/c_boolean"

module Ikra
  Fix        = Struct.new("Fix", :x, :y, :z)
  Volume     = Struct.new("Volume", :x, :y, :z)
  Periodic   = Struct.new("Periodic", :x, :y, :z)
  Coordinate = Struct.new("Coordinate", :x, :y, :z)
  MDLAtom    = Struct.new("MDLAtom", :type, :coordinate, :fix, :visible)
  XYZAtom    = Struct.new("XYZAtom", :type, :coordinate, :energy)

  def summation num
    num == 0 ? 0 : num + summation(num - 1)
  end

  class MDL
    attr_accessor :condition, :volume, :size, :atoms

    def read fp
      @condition     = Condition.new
      @volume, @size = @condition.read_from_mdl(fp)
      @atoms         = []
      fp.each do |line|
        words           = line.split
        atom            = MDLAtom.new
        atom.type       = words[0]
        atom.coordinate = Coordinate.new(Float(words[1]), Float(words[2]), Float(words[3]))
        atom.fix        = Fix.new(Integer(words[4]), Integer(words[5]), Integer(words[6]))
        atom.visible    = Integer(words[7])
        @atoms.push atom
      end
    end

    def write fp
      @condition.write_to_mdl fp, @volume, @size
      @atoms.each do |atom|
        fp.puts "#{atom.type} % 20.15e % 20.15e % 20.15e #{atom.fix.x} #{atom.fix.y} #{atom.fix.z} #{atom.visible}"%[atom.coordinate.x, atom.coordinate.y, atom.coordinate.z]
      end
    end

    def replace_condition c
      # non replaced condition
      c.displacement_ux = Marshal.load(Marshal.dump(@condition.displacement_ux))
      c.displacement_uz = Marshal.load(Marshal.dump(@condition.displacement_uz))
      c.height          = Marshal.load(Marshal.dump(@condition.height))
      c.stress_ex       = Marshal.load(Marshal.dump(@condition.stress_ex))
      c.stress_ez       = Marshal.load(Marshal.dump(@condition.stress_ez))
      c.spbc_dz         = Marshal.load(Marshal.dump(@condition.spbc_dz))

      # deep copy
      @condition = Marshal.load(Marshal.dump(c))
    end
  end

  class XYZ
    attr_accessor :temperature, :atoms, :size
    attr_reader :index_of_frame

    def initialize
      @index_of_frame = -1
    end

    # 1 frame
    def read fp
      @size        = Integer(fp.gets)
      @temperature = Float(fp.gets)
      @atoms       = []
      @size.times do
        words           = fp.gets.split
        atom            = XYZAtom.new
        atom.type       = words[0]
        atom.coordinate = Coordinate.new(Float(words[1]), Float(words[2]), Float(words[3]))
        atom.energy     = Float(words[4])
        @atoms.push atom
      end
      @index_of_frame += 1
    end

    def write fp
      fp.puts @atoms.size, "%20.15e"%[@temperature]
      @atoms.each_with_index do |atom, index|
        fp.puts "#{atom.type} % 20.15e % 20.15e % 20.15e % 20.15e #{index}"%[atom.coordinate.x, atom.coordinate.y, atom.coordinate.z, atom.energy]
      end
    end
  end

  class Element
    attr_accessor :atom_type, :mass
  end

  class Condition
    attr_accessor :loop_number,:quench_steps, :annealing_steps,
      :dynamics_steps, :output_interval, :cg_steps, :dt,
      :temperature, :dtemperature, :elements, :potential_type,
      :displacement, :displacement_u, :displacement_ux, :displacement_uz, :height,
      :shear_stress, :stress, :stress_ex, :stress_ez,
      :spbc_dz, :periodic

    def read fp
      # 全体ループの数
      @loop_number     = Integer(fp.gets)
      # 数値急冷法のステップ数
      @quench_steps    = Integer(fp.gets)
      # 焼き鈍しステップ数（NVTアンサンブル）
      @annealing_steps = Integer(fp.gets)
      # NVEアンサンブルステップ数
      nve              = fp.gets.split
      @dynamics_steps  = Integer(nve[0])
      @output_interval = Integer(nve[1])
      # 共役勾配法ステップ数
      @cg_steps        = Integer(fp.gets)
      # 時間刻みの大きさ（単位：秒）
      @dt              = Float(fp.gets)
      # 設定温度と温度の変化量（焼き鈍しで使う）（単位：K）
      temperatures     = fp.gets.split
      @temperature     = Float(temperatures[0])
      @dtemperature    = Float(temperatures[1])
      # 計算に含まれる原子の種類数
      elements_num     = Integer(fp.gets)
      @elements        = []
      elements_num.times do
        el             = Element.new
        # 原子のラベル（元素記号）と質量（単位：kg）
        words          = fp.gets.split
        el.atom_type   = words[0].chomp
        el.mass        = Float(words[1])
        @elements.push el
      end
      # ポテンシャルタイプの読み込み
      @potential_type  = fp.gets.chomp
      @potential_files = []
      case @potential_type
      when 'pair'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'eam'
        (summation(elements_num) + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'tersoff'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'adp'
        (summation(elements_num) * 3 + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'sw'
        @potential_files.push fp.gets.chomp
      when 'brenner'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when '2bm'
        elements_num.times{ @potential_files.push fp.gets.chomp }
      else
        puts '未定義のポテンシャルタイプです。'
        exit false
      end
      # 変位境界条件（条件を与えるか、その大きさ、与える領域の厚さ）
      disp_bc          = fp.gets.split
      @displacement    = Integer(disp_bc[0]).to_bool
      @displacement_u  = Float(disp_bc[1])
      @displacement_ux = Float(disp_bc[2])
      @displacement_uz = Float(disp_bc[3])
      @height          = Float(disp_bc[4])
      # 応力境界条件（条件を与えるか、その大きさ）
      shear_bc         = fp.gets.split
      @shear_stress    = Integer(shear_bc[0]).to_bool
      @stress          = Float(shear_bc[1])
      @stress_ex       = Float(shear_bc[2])
      @stress_ez       = Float(shear_bc[3])
      # SPBC_dx
      @spbc_dz         = Float(fp.gets)
      # 周期境界条件を与えるか
      periodic_bc      = fp.gets.split
      @periodic        = Periodic.new(
        Integer(periodic_bc[0]).to_bool,
        Integer(periodic_bc[1]).to_bool,
        Integer(periodic_bc[2]).to_bool
      )
    end

    def read_from_mdl fp
      # 全体ループの数
      @loop_number     = Integer(fp.gets)
      # 数値急冷法のステップ数
      @quench_steps    = Integer(fp.gets)
      # 焼き鈍しステップ数（NVTアンサンブル）
      @annealing_steps = Integer(fp.gets)
      # NVEアンサンブルステップ数
      nve = fp.gets.split
      @dynamics_steps  = Integer(nve[0])
      @output_interval = Integer(nve[1])
      # 共役勾配法ステップ数
      @cg_steps        = Integer(fp.gets)
      # 時間刻みの大きさ（単位：秒）
      @dt              = Float(fp.gets)
      # 設定温度と温度の変化量（焼き鈍しで使う）（単位：K）
      temperatures     = fp.gets.split
      @temperature     = Float(temperatures[0])
      @dtemperature    = Float(temperatures[1])
      # 計算に含まれる原子の種類数
      elements_num     = Integer(fp.gets)
      @elements        = []
      elements_num.times do
        el             = Element.new
        # 原子のラベル（元素記号）と質量（単位：kg）
        words          = fp.gets.split
        el.atom_type   = words[0].chomp
        el.mass        = Float(words[1])
        @elements.push el
      end
      # ポテンシャルタイプの読み込み
      @potential_type  = fp.gets.chomp
      @potential_files = []
      case @potential_type
      when 'pair'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'eam'
        (summation(elements_num) + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'tersoff'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'adp'
        (summation(elements_num) * 3 + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'sw'
        @potential_files.push fp.gets.chomp
      when 'brenner'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when '2bm'
        elements_num.times{ @potential_files.push fp.gets.chomp }
      else
        puts '未定義のポテンシャルタイプです。'
        exit false
      end
      # 原子の数
      size             = Integer(fp.gets)
      # 計算領域の大きさ
      volumes          = fp.gets.split
      volume           = Volume.new(
        Float(volumes[0]),
        Float(volumes[1]),
        Float(volumes[2])
      )
      # 変位境界条件（条件を与えるか、その大きさ、与える領域の厚さ）
      disp_bc          = fp.gets.split
      @displacement    = Integer(disp_bc[0]).to_bool
      @displacement_u  = Float(disp_bc[1])
      @displacement_ux = Float(disp_bc[2])
      @displacement_uz = Float(disp_bc[3])
      @height          = Float(disp_bc[4])
      # 応力境界条件（条件を与えるか、その大きさ）
      shear_bc         = fp.gets.split
      @shear_stress    = Integer(shear_bc[0]).to_bool
      @stress          = Float(shear_bc[1])
      @stress_ex       = Float(shear_bc[2])
      @stress_ez       = Float(shear_bc[3])
      # SPBC_dx
      @spbc_dz         = Float(fp.gets)
      # 周期境界条件を与えるか
      periodic_bc      = fp.gets.split
      @periodic        = Periodic.new(
        Integer(periodic_bc[0]).to_bool,
        Integer(periodic_bc[1]).to_bool,
        Integer(periodic_bc[2]).to_bool
      )
      return volume, size
    end

    def write fp
      fp.puts @loop_number,
        @quench_steps,
        @annealing_steps,
        "#{@dynamics_steps} #{@output_interval}",
        @cg_steps,
        "%20.15e"%[@dt],
        "%20.15e %20.15e"%[@temperature, @dtemperature],
        @elements.size
      @elements.each { |el| fp.puts "#{el.atom_type} %20.15e"%[el.mass] }
      fp.puts @potential_type
      @potential_files.each { |pf| fp.puts pf }
      fp.puts "#{@displacement.to_i} %20.15e %20.15e %20.15e %20.15e"%[@displacement_u, @displacement_ux, @displacement_uz, @height],
        "#{@shear_stress.to_i} %20.15e %20.15e %20.15e"%[@stress, @stress_ex, @stress_ez],
        "%20.15e"%[@spbc_dz],
        "#{@periodic.x.to_i} #{@periodic.y.to_i} #{@periodic.z.to_i}"
    end

    def write_to_mdl fp, volume, size
      fp.puts @loop_number,
        @quench_steps,
        @annealing_steps,
        "#{@dynamics_steps} #{@output_interval}",
        @cg_steps,
        "%20.15e"%[@dt],
        "%20.15e %20.15e"%[@temperature, @dtemperature],
        @elements.size
      @elements.each { |el| fp.puts "#{el.atom_type} %20.15e"%[el.mass] }
      fp.puts @potential_type
      @potential_files.each { |pf| fp.puts pf }
      fp.puts size,
        "%20.15e %20.15e %20.15e"%[volume.x, volume.y, volume.z],
        "#{@displacement.to_i} %20.15e %20.15e %20.15e %20.15e"%[@displacement_u, @displacement_ux, @displacement_uz, @height],
        "#{@shear_stress.to_i} %20.15e %20.15e %20.15e"%[@stress, @stress_ex, @stress_ez],
        "%20.15e"%[@spbc_dz],
        "#{@periodic.x.to_i} #{@periodic.y.to_i} #{@periodic.z.to_i}"
    end
  end
end
